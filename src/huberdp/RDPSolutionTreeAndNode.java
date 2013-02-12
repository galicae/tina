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
	
	public String toString() {
		String result = "RDPSolutionTreeAndNode. PA:\n";
		result += this.pa.alignment.getRowAsString(0)+"\n";
		result += this.pa.alignment.getRowAsString(1)+"\n";
		return result;
	}
	
}

/******************************************************************************
 * "A question that sometimes drives me hazy:                                 *
 *  Am I or are the others crazy?"                                            *
 *     - Albert Einstein (1979 - 1966)                                        *
 ******************************************************************************/