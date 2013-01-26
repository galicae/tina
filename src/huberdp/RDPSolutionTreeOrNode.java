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
 * @lastchange 2013-01-26
 *
 */
public class RDPSolutionTreeOrNode
	extends RDPSolutionTreeNode
	implements Comparable<RDPSolutionTreeOrNode> {

	private double score;
		
	public RDPSolutionTreeOrNode() {
		
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
	
}

/******************************************************************************
 * "A question that sometimes drives me hazy:                                 *
 *  Am I or are the others crazy?"                                            *
 *     - Albert Einstein (1979 - 1966)                                        *
 ******************************************************************************/