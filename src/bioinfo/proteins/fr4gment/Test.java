package bioinfo.proteins.fr4gment;

import java.util.HashMap;
import java.util.LinkedList;

import bioinfo.Sequence;
import bioinfo.alignment.SequenceAlignment;
import bioinfo.proteins.PDBEntry;
import bioinfo.proteins.fragm3nt.ProteinFragment;
import bioinfo.superpos.PDBReduce;

public class Test {
	public static void main(String[] args) {
		Sequence seq1 = new Sequence(
				"1nltA01",
				"PQRGKDIKHEISASLEELYKGRTAKLALNKQILVENERKILEVHVEPGMKDGQRIVFKGEADQAPDVIPGDVVFIVSERP");
		Sequence seq2 = new Sequence(
				"1c3gA01",
				"ETVQVNLPVSLEDLFVGKKKSFKIGRKGPHGASEKTQIDIQLKPGWKAGTKITYKNQGDYNPQTGRRKTLQFVIQEKS");

		SequenceAlignment input = new SequenceAlignment(
				seq1,
				seq2,
				"--PQRGKDIKHEIS-----------ASLEELYKGRTAKLALNKQILVENER-----------KILEV-----------------------------HVEPGMKDGQRIVFKGEADQAPDVIPGDVVFIVS----ERP",
				"ET------VQVNLPVSLEDLFVGKK----------------------KSFKIGRKGPHGASEKTQIDIQLKPGWKAGTKITYKNQGDYNPQTGRRK----------------------------TLQFVIQEKS---",
				0.208);

		MultipleSuperposition ms = new MultipleSuperposition(
				"./msp_cluster05/4/cluster88");

		HashMap<Integer, LinkedList<ProteinFragment>> loops = findLoops(ms,
				input);
	}

	/**
	 * this function scans the alignment and returns the start and end points
	 * (in the known sequence, sequence 0) of the aligned regions, thus
	 * identifying the core points
	 * 
	 * @param input
	 *            the alignment (seq0 is the one with known structure and
	 *            multiple superposition)
	 * @return the list of all starts and ends of core segments
	 */
	public static LinkedList<int[]> getCorePoints(SequenceAlignment input) {
		LinkedList<int[]> result = new LinkedList<int[]>();

		int[][] aligned = input.getAlignedResidues();
		int start = 0;
		int end = 0;

		for (int i = 1; i < aligned[0].length; i++) {
			int diff0 = aligned[0][i] - aligned[0][i - 1];
			int diff1 = aligned[1][i] - aligned[1][i - 1];

			if (diff0 == 1 && diff1 == 1) {

			} else {
				end = i - 1;
				int[] temp = { aligned[0][start], aligned[0][end] };
				result.add(temp);
				start = i;
			}
		}
		end = aligned[0].length - 1;
		int[] temp = { aligned[0][start], aligned[0][end] };
		result.add(temp);
		return result;
	}

	/**
	 * this function maps every loop between two core elements on the sequence
	 * of every protein in the superposition and saves the structure fragments
	 * from the superposition corresponding to the loop in a HashMap, the key
	 * being their starting point in the template sequence
	 * 
	 * @param ms
	 *            a multiple superposition where our template sequence/structure
	 *            is included
	 * @param input
	 *            the alignment of the core regions of the query structure and a
	 *            template. Template is at position 0.
	 * @return a HashMap as described above
	 */
	public static HashMap<Integer, LinkedList<ProteinFragment>> findLoops(
			MultipleSuperposition ms, SequenceAlignment input) {
		ms.sort(input.getComponent(0).getID());
		Sequence seq1 = input.getComponent(0);

		double[][] coord = PDBReduce.reduceSinglePDB(ms.getStructures().get(0));
		ProteinFragment usedFrag = new ProteinFragment(seq1.getID(),
				seq1.getSequenceAsString(), coord, coord.length);
		LinkedList<int[]> corePoints = getCorePoints(input);
		HashMap<Integer, LinkedList<ProteinFragment>> betweenCores = new HashMap<Integer, LinkedList<ProteinFragment>>();

		for (int[] loopStart : corePoints) {
			LinkedList<ProteinFragment> currentLoop = new LinkedList<ProteinFragment>();
			for (int i = 0; i < ms.getStructures().size(); i++) {
				PDBEntry pdbX = ms.getStructures().get(i);
				double[][] xCoord = PDBReduce.reduceSinglePDB(pdbX);
				ProteinFragment x = new ProteinFragment(pdbX.getID(),
						pdbX.getSequence(), xCoord, pdbX.length());
				CoreSegmentGotoh got = new CoreSegmentGotoh(-1, -1, 0.5,
						usedFrag, x);

				Sequence xSequence = new Sequence(x.getID(), x.getSequence());

				got.align(seq1, xSequence);
				int[] xCore = got.traceback(corePoints.get(loopStart[0]));
				currentLoop.add(x.getPart(xCore));
			}
			betweenCores.put(loopStart[0], currentLoop);
		}
		return betweenCores;
	}
}
