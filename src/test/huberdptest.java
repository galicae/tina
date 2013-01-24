/******************************************************************************
 * hubeRDPTest is mainly a test routine for calling hubeRDP.                  *
 *                                                                            *
 * This file is best read at line width 80 and tab width 4.                   *
 *                                                                   huberste *
 ******************************************************************************/
package test;

import huberdp.*;

/**
 * @author huberste
 * @lastchange 2013-01-10
 * 
 */
public class huberdptest {

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		RDPSolutionTree t = new RDPSolutionTree();
		RDPPriorityQueue pq = new RDPPriorityQueue(t.getRoot());
		rdp(t, pq);

	}

}

/******************************************************************************
 * "A question that sometimes drives me hazy:                                 *
 *  Am I or are the others crazy?"                                            *
 *     - Albert Einstein (1979 - 1966)                                        *
 ******************************************************************************/