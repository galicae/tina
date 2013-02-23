/******************************************************************************
 * validation.huberdp.HubeRDPCompare.java                                     *
 *                                                                            *
 * Contains the main method to compare against Gotoh.                         *
 *                                                                            *
 * This file is best read at line width 80 and tab width 4.                   *
 *                                                                   huberste *
 ******************************************************************************/
package validation.huberdp;

import java.text.DecimalFormat;
import java.util.LinkedList;
import java.util.Locale;

import bioinfo.Sequence;
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
import huberdp.HubeRDP;
import huberdp.RDPPriorityQueue;
import huberdp.RDPProblem;
import huberdp.RDPSolutionTree;
import huberdp.oracles.*;
import huberdp.scoring.*;
import util.CommandLineParser;

/**
 * HubeRDPCompare compares HubeRDP as an threading tool against Gotoh.
 * 
 * @author huberste
 * @lastchange 2013-02-22
 */
public class HubeRDPCompare {

	private static String USAGE = "Alignment aller Paare in pairfile:\n"
			+ "java HubeRDPCompare --pairs <pairfile> "
			+ "\t--pdbfolder <pdbfilepath>\n\n"
			+ "<pairfile> enthaelt Zeilen der Form \"idone idtwo\"\n\n"
			+ "<pdbfilepath> is the path where the pdbFiles are in."
			+ "Beispiel-Aufruf:\n"
			+ "java -jar HubeRDPValidator.jar --pairs cathscop.inpairs \n"
			+ "\t--pdbfolder ...";

	private static String HEADER = "template\ttarget\trdpscr\trmsd\tgdt\ttmscr\tgotohscr\trmsd\tgdt\ttmscr";

	/**
	 * Validates the HubeRDP.
	 * 
	 * @param args
	 * @throws Exception
	 */
	public static void main(String[] args) {

		// allocate memory
		String pairsString = null, pdbpath = null, hydroFileString = null, ccpFileString = null, voroFileString = null, vpotFileString = null, dsspPathString = null, tempdir = null;
		double gamma = 0.0, delta = 0.0, epsilon = 0.0, zeta = 0.0, gap = 0.0;

		CommandLineParser clp = new CommandLineParser(args);
		if ((pairsString = clp.getStringArg("--pairs")) == null) {
			exitWithUsage("No --pairs was given");
		}
		if ((pdbpath = clp.getStringArg("--pdbfolder")) == null) {
			exitWithUsage("No --pdbfolder was given");
		}
		if ((hydroFileString = clp.getStringArg("--hydro")) == null) {
			exitWithUsage("No --hydro was given");
		}
		if ((ccpFileString = clp.getStringArg("--ccp")) == null) {
			exitWithUsage("No --ccp was given");
		}
		if ((voroFileString = clp.getStringArg("--voro")) == null) {
			exitWithUsage("No --voro was given");
		}
		if ((vpotFileString = clp.getStringArg("--vpot")) == null) {
			exitWithUsage("No --vpot was given");
		}
		if ((dsspPathString = clp.getStringArg("--dssp")) == null) {
			exitWithUsage("No --dssp was given");
		}
		if ((gamma = clp.getDoubleArg("--gamma")) == 0.0) {
			exitWithUsage("No --gamma was given");
		}
		if ((delta = clp.getDoubleArg("--delta")) == 0.0) {
			exitWithUsage("No --delta was given");
		}
		if ((epsilon = clp.getDoubleArg("--epsilon")) == 0.0) {
			exitWithUsage("No --epsilon was given");
		}
		if ((zeta = clp.getDoubleArg("--zeta")) == 0.0) {
			exitWithUsage("no --zeta was given");
		}
		if ((gap = clp.getDoubleArg("--gap")) == 0.0) {
			exitWithUsage("no --gap was given");
		}
		if ((tempdir = clp.getStringArg("--temp")) == null) {
			exitWithUsage("no --temp was given");
		}
		clp = null; // clp -> GC

		// initialize Locale for output
		Locale.setDefault(Locale.US);
		DecimalFormat df = new DecimalFormat("0.0000");

		// initialize PDBFiles
		PDBFileReader fr = new PDBFileReader();

		// initialize HubeRDP
		HubeRDP rdp = new HubeRDP();
		SipplContactPotential sippl = new SipplContactPotential();
		sippl.readFromVPOTFile(vpotFileString);
		HydrophobicityMatrix hydroM = new HydrophobicityMatrix(hydroFileString);
		CCPMatrix ccpM = new CCPMatrix(ccpFileString);
		RDPScoring scoring = new RDPScoring(gamma, delta, epsilon, zeta, gap,
				QuasarMatrix.DAYHOFF_MATRIX, hydroM, ccpM, dsspPathString,
				sippl, null, voroFileString, RDPScoring.GRID_EXTEND,
				RDPScoring.GRID_DENSITY, RDPScoring.GRID_CLASH,
				RDPScoring.MIN_CONTACT, tempdir);
		rdp.setScoring(scoring);
		rdp.addOracle(new RDPOracle(scoring));

		// initialize Gotoh
		LocalSequenceGotoh gotoh = new LocalSequenceGotoh(-10.0, -2.0,
				bioinfo.alignment.matrices.QuasarMatrix.DAYHOFF_MATRIX);

		// initialize TM
		TMMain tmmain = new TMMain();

		// load joblist
		PairFile pairfile = new PairFile(pairsString);
		LinkedList<String[]> joblist = pairfile.getJoblist();

		System.out.println(HEADER);
		
		for (String[] job : joblist) {
			try {
				// load new job
				PDBEntry templateStructure = fr.readPDBFromFile(pdbpath
						+ job[0] + ".pdb");
				Sequence templateSequence = templateStructure.getSequence();
				PDBEntry targetStructure = fr.readPDBFromFile(pdbpath + job[1]
						+ ".pdb");
				Sequence targetSequence = targetStructure.getSequence();

				// align sequences with HubeRDP
				RDPProblem root = new RDPProblem(templateStructure,
						targetSequence);
				RDPSolutionTree t = new RDPSolutionTree(root);
				RDPPriorityQueue pq = new RDPPriorityQueue(t.getRoot());
				rdp.rdp(t, pq);

				// get HubeRDP's (first) alignment
				SequenceAlignment rdpAlignment = t.getRoot().getTA().get(0)
						.getThreading().asSequenceAlignment();

				Transformation rdptr = tmmain.calculateTransformation(
						rdpAlignment, templateStructure, targetStructure);

				SequenceAlignment gotohAlignment = gotoh.align(
						templateSequence, targetSequence);

				Transformation gotohtr = tmmain.calculateTransformation(
						rdpAlignment, templateStructure, targetStructure);

				System.out.println(rdpAlignment.getScore() + "\t"
						+ df.format(rdptr.getRmsd()) + "\t"
						+ df.format(rdptr.getGdt()) + "\t"
						+ df.format(rdptr.getTmscore()) + "\t"
						+ gotohAlignment.getScore() + "\t"
						+ df.format(gotohtr.getRmsd()) + "\t"
						+ df.format(gotohtr.getGdt()) + "\t"
						+ df.format(gotohtr.getTmscore()));
			} catch (Exception e) {
				System.err.println("Error while working on job " + job[0] + " "
						+ job[1]);
				e.printStackTrace();
				System.err.println("continue...");
			}
		}

	}

	public static void exitWithUsage(String error) {
		System.out.println();
		System.out.println(USAGE);
		System.exit(0);
	}

}

/******************************************************************************
 * "A question that sometimes drives me hazy: * Am I or are the others crazy?" *
 * - Albert Einstein (1879 - 1955) *
 ******************************************************************************/
