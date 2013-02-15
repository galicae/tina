/******************************************************************************
 * huberdp.RDPSolutionTreeNode.java                                           *
 * Contains the abstract class RDPSolutionTreeNode which is part of the data  *
 * structure of RDP's solution tree.                                          *
 * It needs to be inherited by AND and OR nodes.                              *
 *                                                                            *
 * This file is best read at line width 80 and tab width 4.                   *
 *                                                                   huberste *
 ******************************************************************************/
package huberdp;

import java.util.LinkedList;

/**
 * RDPSolutionTreeNode is the abstract class that needs to be extended by AND
 * and OR nodes.
 * @author huberste
 * @lastchange 2013-02-11
 */
public abstract class RDPSolutionTreeNode {

	private RDPSolutionTreeNode parent; // for dual link
	private LinkedList<RDPSolutionTreeNode> childs; // for dual link
	
	/**
	 *  TA = TreeAlignment
	 */
	protected LinkedList<TreeAlignment> ta;
	
	/**
	 * states finished if the node can was finished
	 */
	private boolean finished;
	
	/**
	 * costructs a new Node with the given parent
	 * @param parent
	 */
	public RDPSolutionTreeNode(RDPSolutionTreeNode parent) {
		this.setParent(parent);
		this.childs = new LinkedList<RDPSolutionTreeNode>();
		this.setTA(new LinkedList<TreeAlignment>());
		this.setFinished(false);
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
	 * @return the TA
	 */
	public LinkedList<TreeAlignment> getTA() {
		return ta;
	}
	
	/**
	 * @param
	 */
	public void setTA(LinkedList<TreeAlignment> ta) {
		this.ta = ta;
	}
	
	/**
	 * adds a TA to this node.
	 * @param ta
	 */
	public void addTA(TreeAlignment ta) {
		this.ta.add(ta);
	}
	
	/**
	 * Adds multiple TAs to this node.
	 * @param tas
	 */
	public void addTAs(LinkedList<TreeAlignment> tas) {
		for (TreeAlignment ta : tas) {
			this.addTA(ta);
		}
	}
	
	/**
	 * @param finished the finished to set
	 */
	public void setFinished(boolean finished) {
		this.finished = finished;
	}
	
	/**
	 * @return true if node is finished
	 */
	public boolean isFinished() {
		return finished;
	}

	/**
	 * checks if all child are finished
	 * @return true if all child are finished, false else.
	 */
	public boolean checkFinal() {
		// if leaf: final. 
		if (this.isLeaf()) return true;
		
		for (RDPSolutionTreeNode child : childs) {
			if (! child.isFinished()) {
				return false;
			}
		}
		return true;
	}
	
	/**
	 * for RDPSolutionTree.getDepth();
	 * @return the depth of this node's subtree
	 */
	public int getDepth() {
		int result = 0;
		for(RDPSolutionTreeNode child: childs) {
			int temp = child.getDepth();
			if (temp > result)
				result = temp;
//			result = (temp > result) ? temp : result;
		}
		return result+1;
	}
	
}

/******************************************************************************
 * "A question that sometimes drives me hazy:                                 *
 *  Am I or are the others crazy?"                                            *
 *     - Albert Einstein (1879 - 1955)                                        *
 ******************************************************************************/
