/******************************************************************************
 * huberdp.TreeAlignment.java                                                 *
 * Contains the class TreeAlignment which is part of the problem solution of  *
 * RDP.                                                                       *
 *                                                                            *
 * This file is best read at line width 80 and tab width 4.                   *
 *                                                                   huberste *
 ******************************************************************************/
package huberdp;

import bioinfo.Sequence;
import bioinfo.alignment.SequenceAlignment;
import bioinfo.proteins.PDBEntry;

/**
 * TreeAlignment represents the TreeAlignment structure of RDP.
 * @author huberste
 * @lastchange 2013-02-12
 */
public class TreeAlignment extends RDPProblem {
	
	/**
	 * Constructs a new TreeAlignment
	 * @param templateSequence
	 * @param templateStructure
	 * @param targetSequence
	 * @param targetStructure
	 * @param partialAlignment
	 * @param templateStart
	 * @param templateEnd
	 * @param targetStart
	 * @param targetEnd
	 */
	public TreeAlignment(Sequence templateSequence, PDBEntry templateStructure,
			Sequence targetSequence, PDBEntry targetStructure,
			SequenceAlignment partialAlignment,
			int templateStart, int templateEnd, int targetStart, int targetEnd) {
		super(templateSequence, templateStructure,
				targetSequence, targetStructure,
				partialAlignment,
				templateStart, templateEnd, targetStart, targetEnd);
	}
	
	public TreeAlignment(TreeAlignment ta) {
		this(ta.templateSequence, ta.templateStructure,
				ta.targetSequence, ta.targetStructure,
				ta.alignment,
				ta.templateStart, ta.templateEnd, ta.targetStart, ta.targetEnd);
	}
	
	/**
	 * 
	 * @param pa
	 */
	public TreeAlignment(PartialAlignment pa) {
		this(
				pa.templateSequence, pa.templateStructure,
				pa.targetSequence, pa.targetStructure,
				pa.alignment,
				pa.templateStart, pa.templateEnd,
				pa.targetStart, pa.targetEnd);
	}

	/**
	 * toString() method. mostly for debugging.
	 */
	public String toString() {
		String result = "TA\n";
		result += alignment.toStringVerbose();
		return result;
	}
	
	/**
	 * merges two tas.
	 * @param ta1
	 * @param ta2
	 * @return
	 */
	public static TreeAlignment merge(TreeAlignment ta1, TreeAlignment ta2) {
		// TODO
		
		return new TreeAlignment(HubeRDP.mergePaT(ta1, ta2));
/*
		SequenceAlignment alignment = null;
		 
		char[] row10 = ((SequenceAlignment)ta1.alignment).getRow(0);
		char[] row11 = ((SequenceAlignment)ta1.alignment).getRow(1);
		char[] row20 = ((SequenceAlignment)ta2.alignment).getRow(0);
		char[] row21 = ((SequenceAlignment)ta2.alignment).getRow(1);
		
		// break down the alignments into parts (aligned / unaligned)
		
		LinkedList<int[]> parts1 = new LinkedList<int[]>();
		// parts(n)[0] =start,
		// parts(n)[1] = end of part
		LinkedList<Boolean> aligned1 = new LinkedList<Boolean>();
		int n=0; // part

		// begin first part
		parts1.add(new int[2]);
		parts1.get(n)[0] = 0; // set beginning of first part to 0
		aligned1.add(!(row10[n] == '-' || row11[n] == '-')); // set first part
		
		// break down other parts
		for (int pos = 1; pos < row10.length; pos++) {
			if ((row10[pos] == '-' || row11[pos] == '-') == (aligned1.get(n).booleanValue())) {
				parts1.get(n)[1] = pos-1;
				parts1.add(new int[2]);
				aligned1.add(!(row10[pos] == '-' || row11[pos] == '-'));
				n++;
				parts1.get(n)[0] = pos;
			}
		}
		parts1.get(n)[1] = row10.length;
		
		LinkedList<int[]> parts2 = new LinkedList<int[]>();
		// parts(n)[0] =start,
		// parts(n)[1] = end of part
		LinkedList<Boolean> aligned2 = new LinkedList<Boolean>();
		n=0; // part

		// begin first part
		parts2.add(new int[2]);
		parts2.get(n)[0] = 0; // set beginning of first part to 0
		aligned2.add(!(row20[n] == '-' || row21[n] == '-')); // set first part
		
		// break down other parts
		for (int pos = 1; pos < row20.length; pos++) {
			if ((row20[pos] == '-' || row21[pos] == '-') == (aligned2.get(n).booleanValue())) {
				parts2.get(n)[1] = pos-1;
				parts2.add(new int[2]);
				aligned2.add(!(row20[pos] == '-' || row21[pos] == '-'));
				n++;
				parts2.get(n)[0] = pos;
			}
		}
		parts2.get(n)[1] = row20.length;
		n = 0;
		int m = -1;
		// now we have to calculate from position p in row to position m in sequence
		for (int p = 0; p < row10.length; p++) {
			if (row10[p] != '-') {
				m++;
			}
		}
		
		// The following code combines maps. The code seems to be correct, but
		// I can't need it, as I can't recombine an alignment from a map alone.
//		int[][] map2 = ta2.alignment.calcMap();
//		for (int postemp = 0; postemp < map1[0].length; postemp++) {
//			if (map1[0][postemp] == -1) {
//				map1[0][postemp] = map2[0][postemp];
//			}
//		}
//		for (int postar = 0; postar < map1[1].length; postar++) {
//			if (map1[1][postar] == -1) {
//				map1[1][postar] = map2[1][postar];
//			}
//		}
		
		// TODO
		
		String[] rows = new String[2];
		rows[0] += "";
		
		// TODO merge score
		double score = ta1.alignment.getScore();
		

		alignment = new SequenceAlignment(ta1.templateSequence, ta1.targetSequence, rows[0], rows[1], score);
		TreeAlignment result = new TreeAlignment(ta1.templateSequence, ta1.templateStructure, ta1.targetSequence, ta1.targetStructure, alignment);
		return result;
*/
	}
	
}

/******************************************************************************
 * "A question that sometimes drives me hazy:                                 *
 *  Am I or are the others crazy?"                                            *
 *     - Albert Einstein (1879 - 1955)                                        *
 ******************************************************************************/
