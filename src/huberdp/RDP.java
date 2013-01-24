/******************************************************************************
 * rdp is the interface for rdp implementations.                              *
 *                                                                            *
 * This file is best read at line width 80 and tab width 4.                   *
 *                                                                   huberste *
 ******************************************************************************/
package huberdp;

/**
 * @author huberste
 * @lastchange 2013-01-17
 *
 */
public interface RDP {
	
	/**
	 * 
	 * @param t SolutionTree
	 * @param pq PriorityQueue
	 */
	public void rdp(RDPSolutionTree t, RDPPriorityQueue pq);
	
}

/******************************************************************************
 * "A question that sometimes drives me hazy:                                 *
 *  Am I or are the others crazy?"                                            *
 *     - Albert Einstein (1979 - 1966)                                        *
 ******************************************************************************/