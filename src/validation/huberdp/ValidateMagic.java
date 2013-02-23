/******************************************************************************
 * validation.huberdp.ValidateMagic.java                                      *
 *                                                                            *
 * Contains a main method for validating the magic scoring function           *
 *                                                                            *
 * This file is best read at line width 80 and tab width 4.                   *
 *                                                                   huberste *
 ******************************************************************************/
package validation.huberdp;

import files.PairFile;
import huberdp.HubeRDP;
import huberdp.RDPPriorityQueue;
import huberdp.RDPProblem;
import huberdp.RDPSolutionTree;
import huberdp.Scoring;
import huberdp.oracles.RDPOracle;
import huberdp.scoring.MagicScoring;

import java.text.DecimalFormat;
import java.util.LinkedList;
import java.util.Locale;

import util.CommandLineParser;

import bioinfo.alignment.SequenceAlignment;
import bioinfo.alignment.gotoh.LocalSequenceGotoh;
import bioinfo.proteins.PDBEntry;
import bioinfo.proteins.PDBFileReader;
import bioinfo.superpos.TMMain;
import bioinfo.superpos.Transformation;

/**
 * @author huberste
 * 
 */
public class ValidateMagic {

	public static final String USAGE = "Validator for HubeRDP's magicOracle.\n"
			+ "usage:\n"
			+ "java -jar HubeRDPValidator.jar --pairs <pairsfile> --pdb <pdbpath>"
			+ "--tmalign <tmalignpath> --tmp <tmppath>";

	/**
	 * main function for calling HubeRDP
	 * 
	 * @param args
	 *            no args needed
	 * @throws Exception
	 */
	public static void main(String[] args) {

		// allocate memory
		String pairsString = null, pdbpath = null, tmpath, temppath = null;
		// get command Line Arguments
		CommandLineParser clp = new CommandLineParser(args);
		pairsString = clp.getStringArg("--pairs");
		pdbpath = clp.getStringArg("--pdb");
		tmpath = clp.getStringArg("--tmalign");
		temppath = clp.getStringArg("--tmp");
		if (pairsString == null || pdbpath == null || tmpath == null
				|| temppath == null) {
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
		Scoring scoring = new MagicScoring(pdbpath, tmpath, temppath);
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
