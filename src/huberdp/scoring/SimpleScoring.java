/******************************************************************************
 * huberdp.scoring.SimpleScoring.java                                         *
 * Contains the class SimpleScoring that provides a simple scoring function   *
 * for OR nodes.                                                              *
 *                                                                            *
 * This file is best read at line width 80 and tab width 4.                   *
 *                                                                   huberste *
 ******************************************************************************/
package huberdp.scoring;

import huberdp.Scoring;
import bioinfo.alignment.Threading;

/**
 * SimpleScoring provides a simple scoring function for OR nodes.
 * @author huberste
 * @lastchange 2013-02-14
 */
public class SimpleScoring implements Scoring {

	/* (non-Javadoc)
	 * @see huberdp.Scoring#score(huberdp.RDPSolutionTreeOrNode)
	 */
	@Override
	public double score(Threading threading) {
		return (threading.getScore());
	}

	@Override
	public double getScore(Threading t, int m, int n) {
		// TODO Auto-generated method stub
		return 0;
	}

	@Override
	public double getInsertionScore(Threading t, int n) {
		// TODO Auto-generated method stub
		return 0;
	}

	@Override
	public double getDeletionScore(Threading t, int m) {
		// TODO Auto-generated method stub
		return 0;
	}

}

/******************************************************************************
 * "A question that sometimes drives me hazy:                                 *
 *  Am I or are the others crazy?"                                            *
 *     - Albert Einstein (1879 - 1955)                                        *
 ******************************************************************************/
