/**
 * 
 */
package validation.huberdp.sequence;

import files.PairFile;
import huberdp.HubeRDP;
import huberdp.RDPPriorityQueue;
import huberdp.RDPProblem;
import huberdp.RDPSolutionTree;
import huberdp.Scoring;
import huberdp.oracles.RDPOracle;
import huberdp.scoring.SimpleScoring;

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
public class Validator {

	/**
	 * main function for calling HubeRDP
	 * 
	 * @param args
	 *            no args needed
	 * @throws Exception
	 */
	public static void main(String[] args) throws Exception {

		// allocate memory
		String pairsString = null, pdbpath = null;
		// get command Line Arguments
		CommandLineParser clp = new CommandLineParser(args);
		pairsString = clp.getStringArg("--pairs");
		pdbpath = clp.getStringArg("--pdb");
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
		Scoring scoring = new SimpleScoring();
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
				.println("RDP-Scr\tRMSD\tGDT\tTM-Scr\tdepth\tGot-Scr\tRMSD\tGDT\tTM-Scr");

		// INNER LOOP
		for (String[] job : joblist) {
			try {
				// load data
				templateStructure = fr.readPDBFromFile(pdbpath + job[0]
						+ ".pdb");
				targetStructure = fr.readPDBFromFile(pdbpath + job[1]
						+ ".pdb");

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

				System.out.println(df.format(rdpAlignment.getScore()) + "\t"
						+ df.format(rdptmtr.getRmsd()) + "\t"
						+ df.format(rdptmtr.getGdt()) + "\t"
						+ df.format(rdptmtr.getTmscore()) + "\t"
						+ (t.getDepth() / 2) + "\t"
						+ df.format(gotohAlignment.getScore()) + "\t"
						+ df.format(gotohtmtr.getRmsd()) + "\t"
						+ df.format(gotohtmtr.getGdt()) + "\t"
						+ df.format(gotohtmtr.getTmscore()));
			} catch (Exception e) {
				System.err.println("Error occured: " + e.getLocalizedMessage());
				e.printStackTrace();
			}
		}

	}

}
