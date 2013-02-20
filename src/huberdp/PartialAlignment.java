/******************************************************************************
 * huberdp.PartialAlignment.java                                              *
 * Contains the PartialAlignment class which is the definition of an oracle's *
 * partial alignment.                                                                   *
 *                                                                            *
 * This file is best read at line width 80 and tab width 4.                   *
 *                                                                   huberste *
 ******************************************************************************/
package huberdp;

import java.util.LinkedList;

import bioinfo.alignment.SequenceAlignment;
import bioinfo.alignment.Threading;

/**
 * PartialAlignment is an implementation of the Partial Alignment needed for the
 * RDP data structure.
 * 
 * @author huberste
 * @lastchange 2013-02-14
 */
public class PartialAlignment {

	RDPProblem problem;
	SequenceAlignment alignment;

	/**
	 * Constructs a new PartialAlignment given the following parameters:
	 * 
	 * @param problem
	 *            the (sub-)problem this PA solves
	 * @param alignment
	 *            the SequenceAlignment that was given by the oracle
	 */
	public PartialAlignment(RDPProblem problem, SequenceAlignment alignment) {
		this.problem = problem;
		this.alignment = alignment;
	}

	/**
	 * 
	 * @return the SequenceAlignment
	 */
	public SequenceAlignment getAlignment() {
		return alignment;
	}

	/**
	 * 
	 * @return the two subproblems that are defined by the partial alignment
	 */
	public LinkedList<RDPProblem> getSubProblems() {
		LinkedList<RDPProblem> results = new LinkedList<RDPProblem>();

		Threading newThreading = merge();

		char[][] alirows = alignment.getRows();

		// calculate first subproblem
		int pstart = problem.getProblemStart();
		int pend = 0;
		for (int i = 0; i < alignment.length(); i++) {
			if (alirows[0][i] != '-' && alirows[1][i] != '-') {
				pend = problem.getProblemStart() + i;
				break;
			}
		}
		// DONE check sanity of subproblems! (begin < end, ...)
		if (pstart < pend) {
			results.add(new RDPProblem(newThreading, pstart, pend));
		}

		// calculate second subproblem
		pend = problem.getProblemEnd();
		pstart = problem.getProblemEnd();
		for (int i = 0; i < alignment.length(); i++) {
			if (alirows[0][alignment.length() - i - 1] != '-'
					&& alirows[1][alignment.length() - i - i] != '-') {
				pend = problem.getProblemEnd() - i;
				break;
			}
		}
		// DONE check sanity of subproblems! (begin < end, ...)
		if (pstart < pend) {
			results.add(new RDPProblem(newThreading, pstart, pend));
		}

		return results;
	}

	private Threading merge() {
		int[][] oldRows = problem.getThreading().getRows();
		char[][] alignmentRows = alignment.getRows();
		int[][] alipos = alignment.calcPositions();
		int[][] aliRows = new int[2][alignmentRows[0].length];
		// convert from char[][] to int[][]
		int temppos = 0;
		int targpos = 0;
		for (int i = 0; i < alignmentRows[0].length; i++) {
			if (alignmentRows[0][i] == '-') {
				aliRows[0][i] = -1;
			} else {
				aliRows[0][i] = alipos[0][problem.getProblemStart()] + temppos;
				temppos++;
			}
			if (alignmentRows[1][i] == '-') {
				aliRows[1][i] = -1;
			} else {
				aliRows[1][i] = alipos[1][problem.getProblemStart()] + targpos;
				targpos++;
			}
		}

		// calculate newRows[][]
		int[][] newRows = new int[2][];
		newRows[0] = new int[problem.getProblemStart() + alignment.length()
				+ (problem.getThreading().length() - problem.getProblemEnd())];
		newRows[1] = new int[newRows[0].length];
		// copy old start
		for (int i = 0; i < problem.getProblemStart(); i++) {
			newRows[0][i] = oldRows[0][i];
			newRows[1][i] = oldRows[1][i];
		}
		// copy new
		for (int i = 0; i < alignment.length(); i++) {
			newRows[0][problem.getProblemStart() + i] = aliRows[0][i];
			newRows[1][problem.getProblemStart() + i] = aliRows[1][i];
		}
		// copy old end
		for (int i = 0; i < alignment.length(); i++) {
			newRows[0][(problem.getThreading().length() - problem
					.getProblemEnd()) + i] = aliRows[0][(problem.getThreading()
					.length() - problem.getProblemEnd()) + i];
			newRows[1][(problem.getThreading().length() - problem
					.getProblemEnd()) + i] = aliRows[1][(problem.getThreading()
					.length() - problem.getProblemEnd()) + i];
		}

		// FIXME
		double score = 0.0;

		Threading result = new Threading(problem.getThreading().getStructure(),
				problem.getThreading().getSequence(), newRows, score);
		return result;
	}

	/**
	 * toString() function. Mainly for debugging purposes.
	 * 
	 * @return String representation of the partial alignment.
	 */
	public String toString() {
		String result = alignment.getRowAsString(0) + "\n";
		result += alignment.getRowAsString(1);
		return result;
	}

}

/******************************************************************************
 * "A question that sometimes drives me hazy: * Am I or are the others crazy?" *
 * - Albert Einstein (1879 - 1955) *
 ******************************************************************************/
