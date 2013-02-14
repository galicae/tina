/******************************************************************************
 * huberdp.RDPSolutionTree.java                                               *
 * Contains the class RDPSolutionTree which is an implementation of the tree  *
 * data structure of RDP.                                                     *
 *                                                                            *
 * This file is best read at line width 80 and tab width 4.                   *
 *                                                                   huberste *
 ******************************************************************************/
package huberdp;

/**
 * 
 * @author huberste
 * @lastchange 2013-01-17
 */
public class RDPSolutionTree {
	
	private RDPSolutionTreeOrNode root;
	
	/**
	 * constructs an RDPSolutionTree with the root containing the given problem
	 * @param problem
	 */
	public RDPSolutionTree(RDPProblem problem) {
		root = new RDPSolutionTreeOrNode(null, problem, Double.POSITIVE_INFINITY);
	}
	
	/**
	 * 
	 * @return the root of the tree
	 */
	public RDPSolutionTreeOrNode getRoot() {
		return root;
	}
	
}

/******************************************************************************
 * "A question that sometimes drives me hazy:                                 *
 *  Am I or are the others crazy?"                                            *
 *     - Albert Einstein (1979 - 1966)                                        *
 ******************************************************************************/
