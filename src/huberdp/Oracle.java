/******************************************************************************
 * huberdp.Oracle.java                                                        *
 *                                                                            *
 * Contains the abstract class Oracle to be inherited from oracles for a rdp. *
 *                                                                            *
 * This file is best read at line width 80 and tab width 4.                   *
 *                                                                   huberste *
 ******************************************************************************/
package huberdp;

import java.util.LinkedList;

/**
 * Abstract class to be inherited by RDP Oracles.
 * @author huberste
 * @lastchange 2013-02-19
 */
public interface Oracle {
	
	/**
	 * Finds similiar segments, given a certain problem.
	 * @param problem
	 * @param m number of possible solutions the oracle shall find
	 * @return a number of Segments / Alignments
	 */
	public LinkedList<PartialAlignment> findSimiliarSegments
			(RDPProblem problem, int m);
	
}

/******************************************************************************
 * "A question that sometimes drives me hazy:                                 *
 *  Am I or are the others crazy?"                                            *
 *     - Albert Einstein (1879 - 1955)                                        *
 ******************************************************************************/
