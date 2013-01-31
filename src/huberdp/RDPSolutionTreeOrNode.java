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

import java.util.LinkedList;

/**
 * @author huberste
 * @lastchange 2013-01-29
 *
 */
public class RDPSolutionTreeOrNode
	extends RDPSolutionTreeNode
	implements Comparable<RDPSolutionTreeOrNode> {
	
	// SP = SubProblem
	private RDPProblem problem;
	// TA = TreeAlignment
	private LinkedList<TreeAlignment> ta;
	
	private double score;
	
	/**
	 * 
	 * @param parent
	 * @param problem
	 */
	public RDPSolutionTreeOrNode(RDPSolutionTreeAndNode parent, RDPProblem problem) {
		super(parent);
		this.problem=problem;
		this.score = calcScore();
	}

	@Override
	public int compareTo(RDPSolutionTreeOrNode other) {
		if (this.score > other.score) {
			return  1;
		} else if (this.score == other.score) {
			return  0;
		} else { // this.score < other.score
			return -1;
		}
//		return 0;
	}
	
	/**
	 * This function calculates an heuristic score for the Priority Queue
	 * @return
	 */
	private double calcScore() {
		// TODO calculate score
		
		/**
		 * Possible Scoring features:
		 * unclosable gaps
		 * large pags
		 * ...
		 */
		
		return 0.0;
	}

	/**
	 * 
	 * @return the (partial) problem that defines this node
	 */
	public RDPProblem getProblem() {
		return problem;
	}
	
	/**
	 * @return the TreeALignments
	 */
	public LinkedList<TreeAlignment> getTa() {
		return ta;
	}

	/**
	 * @param ta the TreeAlignments to set
	 */
	public void setTa(LinkedList<TreeAlignment> ta) {
		this.ta = ta;
	}
	
	/**
	 * 
	 * @param ta the TreeAlignment to add
	 */
	public void addSolution(TreeAlignment ta) {
		this.ta.add(ta);
	}
	
}

/******************************************************************************
 * "A question that sometimes drives me hazy:                                 *
 *  Am I or are the others crazy?"                                            *
 *     - Albert Einstein (1979 - 1966)                                        *
 ******************************************************************************/