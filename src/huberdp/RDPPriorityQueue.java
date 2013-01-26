/******************************************************************************
 * RDPPriorityQueue is one part of the data structure of RDP                  *
 *                                                                            *
 * This file is best read at line width 80 and tab width 4.                   *
 *                                                                   huberste *
 ******************************************************************************/
package huberdp;

import java.util.PriorityQueue;

/**
 * @author huberste
 * @lastchange 2013-01-17
 *
 */
public class RDPPriorityQueue {

	private PriorityQueue<RDPSolutionTreeOrNode> pq;
	
	/**
	 * Constructs a new (empty) RDPPriorityQueue
	 */
	public RDPPriorityQueue() {
		// TODO correct construction of priorityQueue!
		pq = new PriorityQueue<RDPSolutionTreeOrNode>();
	}
	
	/**
	 * Constructs a new RDPPriorityQueue containing only the given Element
	 * @param node
	 */
	public RDPPriorityQueue(RDPSolutionTreeOrNode node) {
		// TODO correct construction of priorityQueue!
		pq = new PriorityQueue<RDPSolutionTreeOrNode>();
		add(node);
	}
	
	public RDPSolutionTreeOrNode getFirst() {
		return pq.poll();
	}
	
	public boolean isEmpty() {
		return pq.isEmpty();
	}
	
	/**
	 * 
	 * @param e
	 */
	public void add(RDPSolutionTreeOrNode node) {
		pq.add(node);
	}
	
	/**
	 *
	 */
	public void add(RDPSolutionTreeOrNode[] nodes) {
		for (RDPSolutionTreeOrNode node : nodes) {
			pq.add(node);
		}
	}
	
}
