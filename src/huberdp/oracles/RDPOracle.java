/******************************************************************************
 * huberdp.oracles.RDPOracle.java                                             *
 * This file contains the class RDPOracle which is RDP's standard scoring     *
 * function.                                                                  *
 *                                                                            *
 * This file is best read at line width 80 and tab width 4.                   *
 *                                                                   huberste *
 ******************************************************************************/
package huberdp.oracles;

import static bioinfo.alignment.gotoh.Gotoh.FACTOR;
import static util.Util.flip;
import huberdp.Oracle;
import huberdp.PartialAlignment;
import huberdp.RDPProblem;
import huberdp.scoring.RDPScoring;

import java.text.DecimalFormat;
import java.util.LinkedList;

import bioinfo.Sequence;
import bioinfo.alignment.SequenceAlignment;
import bioinfo.alignment.Threading;
import bioinfo.proteins.PDBEntry;

/**
 * @author huberste
 * @lastchange 2013-02-19
 */
public class RDPOracle implements Oracle {

	// for traceback
	protected static final int FROM_TOPLEFT = 1;
	protected static final int FROM_TOP = 2;
	protected static final int FROM_LEFT = 4;

	protected RDPScoring scoring;
	protected RDPProblem problem;
	protected PDBEntry template;
	protected Sequence target;

	public RDPOracle(RDPScoring scoring) {
		this.scoring = scoring;
	}

	/**
	 * @param problem
	 *            the problem we shall find solutions for
	 * @param m
	 *            the maximum number of solutions
	 * @return a LinkedList of local maximum Partial Alignments
	 */
	@Override
	public LinkedList<PartialAlignment> findSimiliarSegments(
			RDPProblem problem, int m) {

		// allocate result
		LinkedList<PartialAlignment> results = new LinkedList<PartialAlignment>();

		this.template = problem.getThreading().getStructure();
		this.target = problem.getThreading().getSequence();

		// set Sequences for Oracle
		String[] rows = problem.getThreading().getRowsAsString();
		String template = "";
		String target = "";
		// only use subproblem: start at ProblemStart, end at ProblemEnd
		for (int i = problem.getProblemStart(); i <= problem.getProblemEnd(); i++) {
			if (rows[0].charAt(i) != '-') {
				template += rows[0].charAt(i);
			}
			if (rows[1].charAt(i) != '-') {
				target += rows[1].charAt(i);
			}
		}
		rows = null; // GC

		// make Sequence Alignment via dynamic programming with scoring
		// function
		SequenceAlignment alignment = align(problem, template, target);

		results.add(new PartialAlignment(problem, alignment));

		return results;
	}

	/**
	 * Simple dynamic programming algorithm.
	 * 
	 * @param templateSequence
	 * @param targetSequence
	 * @return the best local alignment
	 * @TODO make it return several (best) alignments
	 */
	protected SequenceAlignment align(RDPProblem problem, String template,
			String target) {
		// initialize matrices
		int[][] m = new int[template.length() + 1][target.length() + 1];
		int[][] from = new int[template.length() + 1][target.length() + 1];
		Threading threading = problem.getThreading();
		// values for traceback
		int max = Integer.MIN_VALUE;
		int x = 0;
		int y = 0;
		int temppos = problem.getThreading().getFirstAfterInStructure(
				problem.getProblemStart());
		int targpos;
		// fill matrices
		for (int i = 1; i <= template.length(); i++, temppos++) {
			targpos = problem.getThreading().getFirstAfterInSequence(
					problem.getProblemStart());
			for (int j = 1; j <= target.length(); j++, targpos++) {
				// D Matrix and I Matrix are not needed :-)
				// Use D Matrix instead for "comeFrom" values
				int matchValue = m[i - 1][j - 1]
						+ (int) (FACTOR * scoring.getScore(threading,
								temppos, targpos));
				int insertValue = m[i - 1][j]
						+ (int) (FACTOR * scoring.getInsertionScore(threading,
								targpos));
				int deleteValue = m[i][j - 1]
						+ (int) (FACTOR * scoring.getDeletionScore(threading,
								temppos));
				// standard: match
				int mValue = matchValue;
				int fromValue = FROM_TOPLEFT;
				// check Insertion
				if (insertValue == mValue) { // insertion also possible
					fromValue = fromValue | FROM_TOP;
				} else if (insertValue > mValue) { // only insertion
					fromValue = FROM_TOP;
					mValue = insertValue;
				}
				// check deletion
				if (deleteValue == mValue) { // insertion also possible
					fromValue = fromValue | FROM_LEFT;
				} else if (deleteValue > mValue) { // only insertion
					fromValue = FROM_LEFT;
					mValue = deleteValue;
				}
				// local: check 0
				if (mValue < 0) {
					fromValue = 0; // come from no direction
					mValue = 0; // set 0 for local alignment
				}
				// fill M Matrix witch values
				m[i][j] = mValue;
				// fill D Matrix with "comeFrom"
				from[i][j] = fromValue;
				if (mValue >= max) {
					max = mValue;
					x = i - 1;
					y = j - 1;
				}
			}
		}
		
		// begin debugging
//		DecimalFormat df = new DecimalFormat("00.0000");
//		System.out.println("s:\tmax: "+df.format(scoring.smax) + ", min: "+df.format(scoring.smin) + ", avg: " + df.format(scoring.ssum / scoring.num));
//		System.out.println("c:\tmax: "+df.format(scoring.cmax) + ", min: "+df.format(scoring.cmin) + ", avg: " + df.format(scoring.csum / scoring.num));
//		System.out.println("h:\tmax: "+df.format(scoring.hmax) + ", min: "+df.format(scoring.hmin) + ", avg: " + df.format(scoring.hsum / scoring.num));
//		System.out.println("p:\tmax: "+df.format(scoring.pmax) + ", min: "+df.format(scoring.pmin) + ", avg: " + df.format(scoring.psum / scoring.num));
//		System.out.println("result:\tmax: "+df.format(scoring.scoremax) + ", min: "+df.format(scoring.scoremin) + ", avg: " + df.format(scoring.scoresum / scoring.num));
//		System.out.println("insert:\tmax: "+df.format(scoring.insertmax) + ", min: "+df.format(scoring.insertmin) + ", avg: " + df.format(scoring.insertsum / scoring.insertnum));
//		System.out.println("delete:\tmax: "+df.format(scoring.deletemax) + ", min: "+df.format(scoring.deletemin) + ", avg: " + df.format(scoring.deletesum / scoring.deletenum));
//		util.Util.printIntegerArray(m);
//		util.Util.printIntegerArray(from);
		// end debugging

		String row0 = "";
		String row1 = "";
		char actx;
		char acty;

		// start of alignment: unaligned sequences
		for (int i = m[0].length - 1; i > y + 1; i--) {
			row0 += "-";
			row1 += target.charAt(i - 1);
		}
		for (int i = m.length - 1; i > x + 1; i--) {
			row0 += template.charAt(i - 1);
			row1 += "-";
		}

		// real traceback
		while (x >= 0 && y >= 0 && m[x + 1][y + 1] != 0) {
			actx = template.charAt(x);
			acty = target.charAt(y);

			if ((from[x][y] & FROM_TOPLEFT) > 0) { // come from topleft
				row0 += actx;
				row1 += acty;
				y--;
				x--;
			} else if ((from[x][y] & FROM_LEFT) > 0) { // come from left
														// (deletion)
				row0 += "-";
				row1 += acty;
				y--;
			} else if ((from[x][y] & FROM_TOP) > 0) { // come from top
														// (insertion)
				row0 += actx;
				row1 += "-";
				x--;
			}
			else break;
		}
		// end of alignment: unaligned sequences
		for (int i = y + 1; i > 0; i--) {
			row0 += "-";
			row1 += target.charAt(i - 1);
		}
		for (int i = x + 1; i > 0; i--) {
			row0 += template.charAt(i - 1);
			row1 += "-";
		}
		char[] templateRow = flip(row0.toCharArray());
		char[] targetRow = flip(row1.toCharArray());

		SequenceAlignment result = new SequenceAlignment(new Sequence(threading.getStructure()
				.getLongID(), templateRow), new Sequence(threading
				.getSequence().getID(), targetRow), templateRow, targetRow,
				1.0d * max / FACTOR);
		
		// DONE begin debugging
//		System.out.println(result.toStringVerbose());
		// end debugging
		
		return result;

	}
}

/******************************************************************************
 * "A question that sometimes drives me hazy: * Am I or are the others crazy?" *
 * - Albert Einstein (1879 - 1955) *
 ******************************************************************************/
