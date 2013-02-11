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
 * @lastchange 2013-01-29
 *
 */
public abstract class RDPSolutionTreeNode {

	private RDPSolutionTreeNode parent; // for dual link
	private LinkedList<RDPSolutionTreeNode> childs; // for dual link
	
	// TA = TreeAlignment
	protected LinkedList<RDPSolution> ta;
	
	/**
	 * costructs a new Node with the given parent
	 * @param parent
	 */
	public RDPSolutionTreeNode(RDPSolutionTreeNode parent) {
		this.setParent(parent);
		this.childs = new LinkedList<RDPSolutionTreeNode>();
		this.setTA(new LinkedList<RDPSolution>());
	}

	/**
	 * adds a child to the node
	 * @param child
	 */
	public void addChild(RDPSolutionTreeNode child) {
		childs.add(child);
		child.setParent(this);
	}
	
	/**
	 * adds all nodes in the given array to this node's children
	 * @param vSet array containing all the nodes that shall be added to
	 * this node's children
	 */
	public void addChildren(RDPSolutionTreeNode[] vSet) {
		for (RDPSolutionTreeNode node: vSet) {
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
	 * removes the n-th child 
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

	/**
	 * @return the parent
	 */
	public RDPSolutionTreeNode getParent() {
		return parent;
	}

	/**
	 * @param parent the parent to set
	 */
	private void setParent(RDPSolutionTreeNode parent) {
		this.parent = parent;
	}

	/**
	 * @return the ta
	 */
	public LinkedList<RDPSolution> getTA() {
		return ta;
	}
	
	/**
	 * @param
	 */
	public void setTA(LinkedList<RDPSolution> ta) {
		this.ta = ta;
	}
	
	public void addTA(RDPSolution ta) {
		this.ta.add(ta);
	}
	
	public void addTAs(LinkedList<RDPSolution> tas) {
		for (RDPSolution ta : tas) {
			this.addTA(ta);
		}
		
	}
	
}

/******************************************************************************
 * "A question that sometimes drives me hazy:                                 *
 *  Am I or are the others crazy?"                                            *
 *     - Albert Einstein (1979 - 1966)                                        *
 ******************************************************************************/