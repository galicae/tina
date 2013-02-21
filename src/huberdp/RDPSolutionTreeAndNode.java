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

import bioinfo.alignment.Threading;

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

	public void MergeTAs(Scoring scoring) {
		if (isLeaf()) {
			addTA(new TreeAlignment(getPA().getProblem().getThreading()));
		} else {
			LinkedList<RDPSolutionTreeOrNode> children = new LinkedList<RDPSolutionTreeOrNode>();
			for (RDPSolutionTreeNode child : getChilds()) {
				// if (ta.isEmpty()) {
				// addTAs(child.getTA());
				// } else {
				children.add((RDPSolutionTreeOrNode) child);
			}
			LinkedList<TreeAlignment> newTAs = mergeChilds(new LinkedList<TreeAlignment>(), children,
					scoring);
			setTA(newTAs);
		}
	}

	/**
	 * recursive function to merge an AND node's children
	 * 
	 * @return
	 */
	private static LinkedList<TreeAlignment> mergeChilds(
			LinkedList<TreeAlignment> tas,
			LinkedList<RDPSolutionTreeOrNode> children, Scoring scoring) {

		if (children.isEmpty()) {
			return tas;
		}

		LinkedList<TreeAlignment> newTAs;// = new LinkedList<TreeAlignment>();
		LinkedList<RDPSolutionTreeOrNode> remainingChildren = new LinkedList<RDPSolutionTreeOrNode>();

		int problemStart = -1;
		int n = 0;
		// find last SubProblem
		for (int i = 0; i < children.size(); i++) {
			int tempStart = children.get(i).getProblem().getProblemStart();
			if (tempStart > problemStart) {
				problemStart = tempStart;
				n = i;
			}
		}
		// now n is the last subproblem!

		if (tas.isEmpty()) { // No own TAs:
			newTAs = children.get(n).getTA(); // get (last) child's TAs
			for (int i = 0; i < children.size(); i++) {
				if (i != n) {
					remainingChildren.add(children.get(n));
				}
			}
			return mergeChilds(newTAs, remainingChildren, scoring);
		}

		newTAs = new LinkedList<TreeAlignment>();
		RDPSolutionTreeOrNode child = children.get(n);
		RDPProblem problem = child.getProblem();
		LinkedList<TreeAlignment> childsTAs = child.getTA();
		for (TreeAlignment parentTA : tas) {
			for (TreeAlignment childsTA : childsTAs) {

				// Rows of the old (problem)alignment
				int[][] oldRows = parentTA.getThreading().getRows();
				// rows of the new Alignment
				int[][] addRows = childsTA.getThreading().getRows();

				// calculate newRows[][]
				// TODO TEST! THIS CAN'T BE RIGHT I CANT SAY WHY!
				// TODO DOUBLE CHECK!

				int before = problem.getProblemStart();
				int after = parentTA.getThreading().length()
						- problem.getProblemEnd() - 1;
				int afterAtChild = problem.getThreading().length()
						- problem.getProblemEnd() - 1;
				int replace = childsTA.getThreading().length() - before
						- afterAtChild;

				int[][] newRows = new int[2][before + replace + after];
				// copy old start
				for (int i = 0; i < before; i++) {
					newRows[0][i] = oldRows[0][i];
					newRows[1][i] = oldRows[1][i];
				}
				// copy new
				for (int i = 0; i < replace; i++) {
					newRows[0][before + i] = addRows[0][i];
					newRows[1][before + i] = addRows[1][i];
				}
				// copy old end
				for (int i = 0; i < after; i++) {
					newRows[0][before + replace + i] = oldRows[0][problem
							.getProblemEnd() + 1 + i];
					newRows[1][before + replace + i] = oldRows[1][problem
							.getProblemEnd() + 1 + i];
				}

				Threading result = new Threading(parentTA.getThreading()
						.getStructure(), parentTA.getThreading().getSequence(),
						newRows, 0.0);

				// recalculate the score
				double score = scoring.score(result);
				result.setScore(score);
				newTAs.add(new TreeAlignment(result));
			}
		}

		return mergeChilds(newTAs, remainingChildren, scoring);
	}

	/**
	 * toString() method. Mainly for debugging.
	 */
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
