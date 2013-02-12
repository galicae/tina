/******************************************************************************
 * RDPSolutionTreeOrNode is part of the data structure of RDP's solution      *
 * tree.                                                                      *
 *                                                                            *
 * An Or node of a RDPSolutionTree represents a subproblem to be solved.      *
 * An And node of a RDPSolutionTree represents a solution to the supbroblem   *
 * that the nodes father depicts.                                             * 
 *                                                                            *
 * This file is best read at line width 80 and tab width 4.                   *
 *                                                                   huberste *
 ******************************************************************************/
package huberdp;

/**
 * @author huberste
 * @lastchange 2013-02-12
 *
 */
public class RDPSolutionTreeOrNode
	extends RDPSolutionTreeNode
	implements Comparable<RDPSolutionTreeOrNode> {
	
	// SP = SubProblem
	private RDPProblem problem;
	
	private double score;
	
	/**
	 * 
	 * @param parent
	 * @param problem
	 * @param scoring
	 */
	public RDPSolutionTreeOrNode(RDPSolutionTreeAndNode parent, RDPProblem problem, Scoring scoring) {
		super(parent);
		this.problem=problem;
		this.setScore(scoring.score(this));
	}
	
	/**
	 * 
	 * @param parent
	 * @param problem
	 * @param score
	 */
	public RDPSolutionTreeOrNode(RDPSolutionTreeAndNode parent, RDPProblem problem, double score) {
		super(parent);
		this.problem=problem;
		this.setScore(score);
	}

	@Override
	public int compareTo(RDPSolutionTreeOrNode other) {
		if (this.getScore() > other.getScore()) {
			return  1;
		} else if (this.getScore() == other.getScore()) {
			return  0;
		} else { // this.score < other.score
			return -1;
		}
//		return 0;
	}

	/**
	 * 
	 * @return the (sub-) problem that defines this node
	 */
	public RDPProblem getProblem() {
		return problem;
	}
	
	/**
	 * simple toString() function, mainly for debugging
	 * @return String representation of this OR node
	 */
	public String toString() {
		String result = "RDPSolutionTreeOrNode. ";

		if (this.problem.alignment != null) {
			result += "Partial alignment up to now:\n"+
					this.problem.alignment.toStringVerbose()+
					"\n";
		} 
		result += "SP:\n" +this.problem.toString();
		return result;
	}

	/**
	 * @return the score
	 */
	public double getScore() {
		return score;
	}

	/**
	 * @param score the score to set
	 */
	private void setScore(double score) {
		this.score = score;
	}
	
}

/******************************************************************************
 * "A question that sometimes drives me hazy:                                 *
 *  Am I or are the others crazy?"                                            *
 *     - Albert Einstein (1979 - 1966)                                        *
 ******************************************************************************/