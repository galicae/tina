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

	private RDPProblem problem;
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
		this.problem =problem;
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
	public LinkedList<RDPProblem> getSubProblems(Scoring scoring) {
		LinkedList<RDPProblem> results = new LinkedList<RDPProblem>();

		Threading newThreading = merge(scoring);

		char[][] alirows = alignment.getRows();

		// calculate first subproblem
		int pstart = getProblem().getProblemStart();
		int pend = 0;
		for (int i = 0; i < alignment.length(); i++) {
			if (alirows[0][i] != '-' && alirows[1][i] != '-') {
				pend = getProblem().getProblemStart() + i-1;
				break;
			}
		}
		// DONE check sanity of subproblems! (begin < end, ...)
		if (pstart < pend) {
			results.add(new RDPProblem(newThreading, pstart, pend));
		}

		// calculate second subproblem
		pend = getProblem().getProblemEnd();
		pstart = getProblem().getProblemEnd();
		for (int i = 0; i < alignment.length(); i++) {
			if (alirows[0][alignment.length() - i - 1] != '-'
					&& alirows[1][alignment.length() - i - 1] != '-') {
				pend = getProblem().getProblemEnd() - i;
				break;
			}
		}
		// DONE check sanity of subproblems! (begin < end, ...)
		if (pstart < pend) {
			results.add(new RDPProblem(newThreading, pstart, pend));
		}

		return results;
	}

	private Threading merge(Scoring scoring) {
		// Rows of the old (problem)alignment
		int[][] oldRows = getProblem().getThreading().getRows();
		// rows of the new Alignment
		int[][] aliRows = alignment.getRowsAsIntArray();

		// calculate newRows[][]
		int[][] newRows = new int[2][];
		newRows[0] = new int[getProblem().getProblemStart() + alignment.length()
				+ (getProblem().getThreading().length() - getProblem().getProblemEnd())
				- 1];
		newRows[1] = new int[newRows[0].length];
		// copy old start
		for (int i = 0; i < getProblem().getProblemStart(); i++) {
			newRows[0][i] = oldRows[0][i];
			newRows[1][i] = oldRows[1][i];
		}
		// copy new
		for (int i = getProblem().getProblemStart(); i < getProblem().getProblemStart() + alignment.length(); i++) {
			newRows[0][i] = aliRows[0][i];
			newRows[1][i] = aliRows[1][i];
		}
		// copy old end
		for (int i = getProblem().getProblemStart() + alignment.length(); i < newRows[0].length; i++) {
			newRows[0][i] = oldRows[0][(getProblem().getProblemEnd()) - getProblem().getProblemStart() - alignment.length() + i + 1];
			newRows[1][i] = oldRows[1][(getProblem().getProblemEnd()) - getProblem().getProblemStart() - alignment.length() + i + 1];
		}

		Threading result = new Threading(getProblem().getThreading().getStructure(),
				getProblem().getThreading().getSequence(), newRows, 0.0);

		// recalculate the score
		double score = scoring.score(result);
		result.setScore(score);

		
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

	/**
	 * @return the problem
	 */
	public RDPProblem getProblem() {
		return problem;
	}

}

/******************************************************************************
 * "A question that sometimes drives me hazy: * Am I or are the others crazy?" *
 * - Albert Einstein (1879 - 1955) *
 ******************************************************************************/
