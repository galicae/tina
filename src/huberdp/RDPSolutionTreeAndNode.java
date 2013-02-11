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
 * @lastchange 2013-02-11
 *
 */
public class RDPSolutionTreeAndNode extends RDPSolutionTreeNode {

	// PA = PartialAlignment
	private PartialAlignment pa;
	
	/**
	 * Constructs an RDPSolutionTreeAndNode
	 * @param parent the node's parent
	 * @param alignment the Partial alignment that defines this node
	 */
	public RDPSolutionTreeAndNode(RDPSolutionTreeOrNode parent, PartialAlignment alignment) {
		super(parent);
		this.setAlignment(alignment);
	}

	/**
	 * @return the alignment
	 */
	public PartialAlignment getPA() {
		return pa;
	}

	/**
	 * @param pa the alignment to set
	 */
	public void setAlignment(PartialAlignment pa) {
		this.pa = pa;
	}
	
	/**
	 * merges the given TA with this AND node's PA and adds it to this node's
	 * list of TAs
	 * @param ta the ta to merge
	 */
	public void mergeTA(RDPSolution ta) {
		// TODO merge ta with this.alignment
		this.addTA(ta);
	}
	
}

/******************************************************************************
 * "A question that sometimes drives me hazy:                                 *
 *  Am I or are the others crazy?"                                            *
 *     - Albert Einstein (1979 - 1966)                                        *
 ******************************************************************************/