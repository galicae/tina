/******************************************************************************
 * huberdp.RDPSolutionTreeOrNode.java                                         *
 *                                                                            *
 * Contains the class RDPSolutionTreeOrNode which is part of the data         *
 * structure of RDP's solution tree.                                          *
 *                                                                            *
 * An OR node of a RDPSolutionTree represents a (sub-)problem to be solved.   *
 * An AND node of a RDPSolutionTree represents a solution to the supbroblem   *
 * that the nodes father depicts.                                             * 
 *                                                                            *
 * This file is best read at line width 80 and tab width 4.                   *
 *                                                                   huberste *
 ******************************************************************************/
package huberdp;

/**
 * RDPSolutionTreeOrNode represents a (sub-)problem to be solved.
 * 
 * @author huberste
 * @lastchange 2013-02-12
 */
public class RDPSolutionTreeOrNode extends RDPSolutionTreeNode implements
		Comparable<RDPSolutionTreeOrNode> {

	// SP = SubProblem
	private RDPProblem problem;

	private double score;
	
	/**
	 * Standard Constructor
	 * @param parent
	 * @param problem
	 * @param score
	 */
	public RDPSolutionTreeOrNode(RDPSolutionTreeAndNode parent,
			RDPProblem problem, double score) {
		super(parent);
		this.problem = problem;
		this.setScore(score);
	}
	
	/**
	 * 
	 * @param parent
	 * @param problem
	 * @param scoring
	 */
	public RDPSolutionTreeOrNode(RDPSolutionTreeAndNode parent,
			RDPProblem problem, Scoring scoring) {
		this(parent, problem, scoring.score(problem.getThreading()));
	}
	
	@Override
	public int compareTo(RDPSolutionTreeOrNode other) {
		if (this.getScore() > other.getScore()) {
			return 1;
		} else if (this.getScore() == other.getScore()) {
			return 0;
		} else { // this.score < other.score
			return -1;
		}
		// return 0;
	}

	/**
	 * 
	 * @return the (sub-) problem that defines this node
	 */
	public RDPProblem getProblem() {
		return problem;
	}

	/**
	 * @return the score
	 */
	public double getScore() {
		return score;
	}

	/**
	 * @param score
	 *            the score to set
	 */
	private void setScore(double score) {
		this.score = score;
	}

	/**
	 * simple toString() function, mainly for debugging
	 * 
	 * @return String representation of this OR node
	 */
	public String toString() {
		String result = "RDPSolutionTreeOrNode. ";

		if (this.problem.getThreading() != null) {
			result += "Partial alignment up to now:\n"
					+ this.problem.getThreading().toStringVerbose() + "\n";
		}
		result += "SP:\n" + this.problem.toString();
		return result;
	}

}

/******************************************************************************
 * "A question that sometimes drives me hazy:                                 *
 *  Am I or are the others crazy?"                                            *
 * - Albert Einstein (1879 - 1955)                                            *
 ******************************************************************************/
