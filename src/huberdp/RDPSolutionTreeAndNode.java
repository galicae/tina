/******************************************************************************
 * RDPSolutionTreeOrNode is part of the data structure of RDP's solution      *
 * tree.                                                                      *
 *                                                                            *
 * An And node of a RDPSolutionTree represents a solution to the supbroblem   *
 * that the nodes father depicts.                                             * 
 * An Or node of a RDPSolutionTree represents a subproblem to be solved.      *
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
public class RDPSolutionTreeAndNode extends RDPSolutionTreeNode {

	// PA = PartialAlignment
	private RDPProblem alignment;
	// TA = TreeAlignment
	private TreeAlignment ta;
	
	/**
	 * Constructs an RDPSolutionTreeAndNode
	 * @param parent the node's parent
	 * @param alignment the Partial alignment that defines this node
	 */
	public RDPSolutionTreeAndNode(RDPSolutionTreeOrNode parent, RDPProblem alignment) {
		super(parent);
		this.setAlignment(alignment);
	}

	/**
	 * @return the alignment
	 */
	public RDPProblem getAlignment() {
		return alignment;
	}

	/**
	 * @param alignment the alignment to set
	 */
	public void setAlignment(RDPProblem alignment) {
		this.alignment = alignment;
	}

	/**
	 * @return the TreeAlignment
	 */
	public TreeAlignment getTa() {
		return ta;
	}

	/**
	 * @param ta the TreeAlignment to set
	 */
	public void setTa(TreeAlignment ta) {
		this.ta = ta;
	}
	
}

/******************************************************************************
 * "A question that sometimes drives me hazy:                                 *
 *  Am I or are the others crazy?"                                            *
 *     - Albert Einstein (1979 - 1966)                                        *
 ******************************************************************************/