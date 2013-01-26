/******************************************************************************
 * RDPSolutionTreeNode is part of the data structure of RDP's solution tree.  *
 *                                                                            *
 * This file is best read at line width 80 and tab width 4.                   *
 *                                                                   huberste *
 ******************************************************************************/
package huberdp;

import java.util.LinkedList;

/**
 * @author huberste
 *
 */
public abstract class RDPSolutionTreeNode {

	private LinkedList<RDPSolutionTreeNode> childs;
	
	/**
	 * adds a child to the node
	 * @param child
	 */
	public void addChild(RDPSolutionTreeNode child) {
		childs.add(child);
	}
	
	/**
	 * adds all nodes in the given array to this node's children
	 * @param children array containing all the nodes that shall be added to
	 * this node's children
	 */
	public void addChildren(RDPSolutionTreeNode[] children) {
		for (RDPSolutionTreeNode node: children) {
			this.addChild(node);
		}
	}
	
	/**
	 * 
	 * @param arg a RDPSolutionTreeNode
	 * @return true if arg is a child of this node, false else
	 */
	public boolean contains(RDPSolutionTreeNode arg) {
		for (RDPSolutionTreeNode child : childs) {
			if (child.equals(arg)) return true;
		}
		return false;
	}
	
	/**
	 * 
	 */
	public void removeChild(int n) {
		childs.remove(n);
	}
	
	/**
	 * @return the childs
	 */
	public LinkedList<RDPSolutionTreeNode> getChilds() {
		return childs;
	}

	/**
	 * 
	 * @return true if this node has no childs
	 */
	public boolean isLeaf() {
		if (childs.size() == 0) {
			return true;
		} else {
			return false;
		}
	}
	
}

/******************************************************************************
 * "A question that sometimes drives me hazy:                                 *
 *  Am I or are the others crazy?"                                            *
 *     - Albert Einstein (1979 - 1966)                                        *
 ******************************************************************************/