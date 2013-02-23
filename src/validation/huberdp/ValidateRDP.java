/******************************************************************************
 * validation.huberdp.ValidateRDP.java                                        *
 *                                                                            *
 * Contains a main method for validating the RDPScoring function.             *
 *                                                                            *
 * This file is best read at line width 80 and tab width 4.                   *
 *                                                                   huberste *
 ******************************************************************************/
package validation.huberdp;

import huberdp.HubeRDP;
import huberdp.RDPPriorityQueue;
import huberdp.RDPProblem;
import huberdp.RDPSolutionTree;
import huberdp.Scoring;
import huberdp.oracles.RDPOracle;
import huberdp.scoring.RDPScoring;

import java.text.DecimalFormat;
import java.util.LinkedList;
import java.util.Locale;

import util.CommandLineParser;
import bioinfo.alignment.SequenceAlignment;
import bioinfo.alignment.gotoh.LocalSequenceGotoh;
import bioinfo.alignment.matrices.QuasarMatrix;
import bioinfo.energy.potential.SipplContactPotential;
import bioinfo.energy.potential.hydrophobicity.HydrophobicityMatrix;
import bioinfo.proteins.CCPMatrix;
import bioinfo.proteins.PDBEntry;
import bioinfo.proteins.PDBFileReader;
import bioinfo.superpos.TMMain;
import bioinfo.superpos.Transformation;
import files.PairFile;

/**
 * @author huberste
 * @lastchange 2013-02-23
 */
public class ValidateRDP {

	public static final String USAGE = "Validator for HubeRDP's magicOracle.\n"
			+ "usage:\n"
			+ "java -jar HubeRDPValidator.jar --pairs <pairsfile> "
			+ "--pdb <pdbpath> --tmp <temppath> "
			+ "--hydro <hydrofile> --ccp <ccpfile> --voro <voro++> "
			+ "--vpot <vpotfile> --dssp <dssppath> "
			+ "--gamma <gamma> --delta <delta> --epsilon <epsilon> "
			+ "--zeta <zeta> --gap <gap>";

	/**
	 * main function for calling HubeRDP
	 * 
	 * @param args
	 *            no args needed
	 * @throws Exception
	 */
	public static void main(String[] args) {

		// allocate memory
		String pairsString = null, pdbpath = null, temppath = null, hydroFileString = null, ccpFileString = null, voroFileString = null, vpotFileString = null, dsspPathString = null;
		double gamma = 0.0, delta = 0.0, epsilon = 0.0, zeta = 0.0, gap = 0.0;
		// get command Line Arguments
		CommandLineParser clp = new CommandLineParser(args);
		pairsString = clp.getStringArg("--pairs");
		pdbpath = clp.getStringArg("--pdb");
		temppath = clp.getStringArg("--tmp");
		hydroFileString = clp.getStringArg("--hydro");
		ccpFileString = clp.getStringArg("--ccp");
		voroFileString = clp.getStringArg("--voro");
		vpotFileString = clp.getStringArg("--vpot");
		dsspPathString = clp.getStringArg("--dssp");
		gamma = clp.getDoubleArg("--gamma");
		delta = clp.getDoubleArg("--delta");
		epsilon = clp.getDoubleArg("--epsilon");
		zeta = clp.getDoubleArg("--zeta");
		gap = clp.getDoubleArg("--gap");

		if (pairsString == null || pdbpath == null || temppath == null
				|| hydroFileString == null || ccpFileString == null
				|| voroFileString == null || vpotFileString == null
				|| dsspPathString == null || gamma == 0.0 || delta == 0.0
				|| epsilon == 0.0 || zeta == 0.0 || gap == 0.0) {
			System.out.println(USAGE);
			System.exit(0);
		}
		clp = null; // clp -> GC

		// load joblist
		PairFile pairfile = new PairFile(pairsString);
		LinkedList<String[]> joblist = pairfile.getJoblist();

		// set test data
		PDBEntry templateStructure = null;
		PDBEntry targetStructure = null;

		// initialize PDBFileReader
		PDBFileReader fr = new PDBFileReader();

		// construct RDP
		HubeRDP rdp = new HubeRDP();
		// set scoring
		SipplContactPotential sippl = new SipplContactPotential();
		sippl.readFromVPOTFile(vpotFileString);
		Scoring scoring = new RDPScoring(gamma, delta, epsilon, zeta, gap,
				QuasarMatrix.DAYHOFF_MATRIX, new HydrophobicityMatrix(
						hydroFileString), new CCPMatrix(ccpFileString),
				dsspPathString, sippl, templateStructure, voroFileString,
				RDPScoring.GRID_EXTEND, RDPScoring.GRID_DENSITY,
				RDPScoring.GRID_CLASH, RDPScoring.MIN_CONTACT, temppath);
		rdp.setScoring(scoring);
		// add oracles
		rdp.addOracle(new RDPOracle(scoring));

		LocalSequenceGotoh gotoh = new LocalSequenceGotoh(-10.0, -2.0,
				bioinfo.alignment.matrices.QuasarMatrix.DAYHOFF_MATRIX);

		// initialize TM stuff
		TMMain tmmain = new TMMain();

		// initialize output stuff
		Locale.setDefault(Locale.US);
		DecimalFormat df = new DecimalFormat("0.0000");

		System.out
				.println("tmplt\ttrgt\tRDP-Scr\tRMSD\tGDT\tTM-Scr\tdpth\tGth-Scr\tRMSD\tGDT\tTM-Scr");

		// INNER LOOP
		for (String[] job : joblist) {
			try {
				// load data
				templateStructure = fr.readPDBFromFile(pdbpath + job[0]
						+ ".pdb");
				targetStructure = fr.readPDBFromFile(pdbpath + job[1] + ".pdb");

				// construct rdp tree
				RDPProblem root = new RDPProblem(templateStructure,
						targetStructure.getSequence());
				RDPSolutionTree t = new RDPSolutionTree(root);
				// construct priority queue
				RDPPriorityQueue pq = new RDPPriorityQueue(t.getRoot());
				// execute rdp algorithm
				rdp.rdp(t, pq);
				// Solutions are now in t.getRoot();

				// get HubeRDP's (first) alignment
				SequenceAlignment rdpAlignment = t.getRoot().getTA().get(0)
						.getThreading().asSequenceAlignment();

				SequenceAlignment gotohAlignment = gotoh.align(
						templateStructure.getSequence(),
						targetStructure.getSequence());

				Transformation rdptmtr = tmmain.calculateTransformation(
						rdpAlignment, templateStructure, targetStructure);

				Transformation gotohtmtr = tmmain.calculateTransformation(
						gotohAlignment, templateStructure, targetStructure);

				System.out.println(job[0] + "\t" + job[1] + "\t"
						+ df.format(rdpAlignment.getScore()) + "\t"
						+ df.format(rdptmtr.getRmsd()) + "\t"
						+ df.format(rdptmtr.getGdt()) + "\t"
						+ df.format(rdptmtr.getTmscore()) + "\t"
						+ (t.getDepth() / 2) + "\t"
						+ df.format(gotohAlignment.getScore()) + "\t"
						+ df.format(gotohtmtr.getRmsd()) + "\t"
						+ df.format(gotohtmtr.getGdt()) + "\t"
						+ df.format(gotohtmtr.getTmscore()));
			} catch (Exception e) {
				System.err.println("Error occured at " + job[0] + " " + job[1]
						+ ": " + e.getLocalizedMessage());
				e.printStackTrace();
				System.err.println("continue...");
			}
		}

	}

}

/******************************************************************************
 * "A question that sometimes drives me hazy: * Am I or are the others crazy?" *
 * - Albert Einstein (1879 - 1955) *
 ******************************************************************************/
