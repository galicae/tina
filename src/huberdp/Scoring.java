/******************************************************************************
 * Scoring is an abstract class for scoring OR-Nodes                          *
 *                                                                            *
 * This file is best read at line width 80 and tab width 4.                   *
 *                                                                   huberste *
 ******************************************************************************/
package huberdp;

/**
 * Scoring is an abstract class for scoring OR-Nodes
 *
 * Possible Scoring features:
 * * unclosable gaps
 * * large gaps
 * ...
 * 
 * @author huberste
 * @lastchange 2013-02-12
 */
public interface Scoring {

	/**
	 * scores an OR node
	 * @param node
	 * @return
	 */
	public double score(RDPSolutionTreeOrNode node);
	
}

/******************************************************************************
 * "A question that sometimes drives me hazy:                                 *
 *  Am I or are the others crazy?"                                            *
 *     - Albert Einstein (1979 - 1966)                                        *
 ******************************************************************************/