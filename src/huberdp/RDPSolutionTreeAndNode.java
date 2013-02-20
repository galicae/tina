/******************************************************************************
 * huberdp.RDPSolutionTreeAndNode.java                                        *
 * Contains the class RDPSolutionTreeAndNode which is part of the tree data   *
 * structure of RDP's solution tree.                                          *
 *                                                                            *
 * An AND node of a RDPSolutionTree represents a solution to the supbroblem   *
 * that the node's father depicts.                                            * 
 * An OR node of a RDPSolutionTree represents a (sub-)problem to be solved.   *
 *                                                                            *
 * This file is best read at line width 80 and tab width 4.                   *
 *                                                                   huberste *
 ******************************************************************************/
package huberdp;

import java.util.LinkedList;

/**
 * RDPSolutionTreeAndNode represents a solution to the supbroblem that the
 * node's father depicts.
 * 
 * @author huberste
 * @lastchange 2013-02-12
 */
public class RDPSolutionTreeAndNode extends RDPSolutionTreeNode {

	// PA = PartialAlignment
	private PartialAlignment pa;

	/**
	 * Constructs an RDPSolutionTreeAndNode
	 * 
	 * @param parent
	 *            the node's parent
	 * @param alignment
	 *            the Partial alignment that defines this node
	 */
	public RDPSolutionTreeAndNode(RDPSolutionTreeOrNode parent,
			PartialAlignment alignment) {
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
	 * @param pa
	 *            the alignment to set
	 */
	public void setAlignment(PartialAlignment pa) {
		this.pa = pa;
	}

	public void MergeTAs() {
		if (isLeaf()) {
			addTA(new TreeAlignment(getPA().getProblem().getThreading()));
		} else {
			for (RDPSolutionTreeNode child : getChilds()) {
				if (ta.isEmpty()) {
					addTAs(child.getTA());
				} else {
					LinkedList<TreeAlignment> newTAs = new LinkedList<TreeAlignment>();
					for (TreeAlignment ta1 : child.getTA()) {
						for (TreeAlignment ta2 : getTA()) {
							// TODO MERGE TAs
							// oldsource:
							// newTAs.add(TreeAlignment.merge(ta1,
							// ta2));
						}
					}
					setTA(newTAs);
				}
			}
		}
	}

	public String toString() {
		String result = "RDPSolutionTreeAndNode. PA:\n";
		result += this.pa.alignment.getRowAsString(0) + "\n";
		result += this.pa.alignment.getRowAsString(1) + "\n";
		return result;
	}

}

/******************************************************************************
 * "A question that sometimes drives me hazy: * Am I or are the others crazy?" *
 * - Albert Einstein (1879 - 1955) *
 ******************************************************************************/
