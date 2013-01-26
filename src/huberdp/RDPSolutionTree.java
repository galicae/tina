/******************************************************************************
 * RDPSolutionTree is one part of the the data structure of RDP.              *
 *                                                                            *
 * This file is best read at line width 80 and tab width 4.                   *
 *                                                                   huberste *
 ******************************************************************************/
package huberdp;

/**
 * @author huberste
 * @lastchange 2013-01-17
 *
 */
public class RDPSolutionTree extends RDPSolutionTreeOrNode {
	
	private RDPSolutionTreeOrNode root;
	
	/**
	 * Constructs an empty RDPSolutionTree
	 */
	public RDPSolutionTree() {
		root = null;
	}
	
	/**
	 * constructs an RDPSolutionTree with the root containing the given problem
	 * @param problem
	 */
	public RDPSolutionTree(RDPProblem problem) {
		root = new RDPSolutionTreeOrNode(problem);
		
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