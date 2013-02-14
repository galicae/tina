/******************************************************************************
 * huberdp.RDP                                                                *
 * Contains the interface RDP for RDP implementations.                        *
 *                                                                            *
 * This file is best read at line width 80 and tab width 4.                   *
 *                                                                   huberste *
 ******************************************************************************/
package huberdp;

/**
 * RDP is an interface to be used by every implementation of RDP.
 * @author huberste
 * @lastchange 2013-02-14
 */
public interface RDP {
	
	/**
	 * Construct a new RDP
	 * @param t SolutionTree
	 * @param pq PriorityQueue
	 */
	public RDPSolutionTreeNode rdp(RDPSolutionTree t, RDPPriorityQueue pq);
	
}

/******************************************************************************
 * "A question that sometimes drives me hazy:                                 *
 *  Am I or are the others crazy?"                                            *
 *     - Albert Einstein (1879 - 1955)                                        *
 ******************************************************************************/
