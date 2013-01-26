/******************************************************************************
 * Oracle is the interface for oracles for a rdp.                             *
 *                                                                            *
 * This file is best read at line width 80 and tab width 4.                   *
 *                                                                   huberste *
 ******************************************************************************/
package huberdp;


import java.util.LinkedList;

/**
 * @author huberste
 * @lastchange 2013-01-17
 */
public interface Oracle {
	
	/**
	 * 
	 * @param problem
	 * @param m number of possible solutions the oracle shall find
	 * @return a number of Segments / Alignments
	 */
	public LinkedList<RDPProblem> findSimiliarSegments(RDPProblem problem, int m);
	
}

/******************************************************************************
 * "A question that sometimes drives me hazy:                                 *
 *  Am I or are the others crazy?"                                            *
 *     - Albert Einstein (1979 - 1966)                                        *
 ******************************************************************************/