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

	private PriorityQueue<RDPSolutionTree> pq;
	
	public RDPPriorityQueue(RDPSolutionTreeNode root) {
		// TODO correct construction of priorityQueue!
		pq = new PriorityQueue<RDPSolutionTree>();
	}
	
	public RDPSolutionTree getFirst() {
		return pq.poll();
	}
	
}
