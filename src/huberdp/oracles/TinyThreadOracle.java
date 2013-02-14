/******************************************************************************
 * huberdp.Oracles.TinyThreadOracle.java                                      *
 * This file contains the class TinyThreadOracle which is an oracle for rdp   *
 * making threading possible.                                                 *
 *                                                                            *
 * This file is best read at line width 80 and tab width 4.                   *
 *                                                                   huberste *
 ******************************************************************************/
package huberdp.oracles;

import java.util.LinkedList;

import huberdp.Oracle;
import huberdp.PartialAlignment;
import huberdp.RDPProblem;

/**
 * TinyThreadOracle is an oracle for HubeRDP that makes threading possible.
 * @author huberste
 * @lastchange 2013-02-14
 */
public class TinyThreadOracle extends Oracle {

	/**
	 * the function that is called by gAND.
	 * It is the main function of the oracle.
	 * @author huberste
	 * @param problem a (sub-)problem for which PartialAlignments are to be found
	 * @param m maximum number of PartialAlignments to be returned
	 * @return a LinkedList of at most m PartialAlignments
	 */
	@Override
	public LinkedList<PartialAlignment> findSimiliarSegments
	(RDPProblem problem, int m)	{
		return null;
	}

}

/******************************************************************************
 * "A question that sometimes drives me hazy:                                 *
 *  Am I or are the others crazy?"                                            *
 *     - Albert Einstein (1879 - 1955)                                        *
 ******************************************************************************/
