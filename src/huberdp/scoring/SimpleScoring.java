/******************************************************************************
 * huberdp.scoring.SimpleScoring.java                                         *
 * Contains the class SimpleScoring that provides a simple scoring function   *
 * for OR nodes.                                                              *
 *                                                                            *
 * This file is best read at line width 80 and tab width 4.                   *
 *                                                                   huberste *
 ******************************************************************************/
package huberdp.scoring;

import bioinfo.alignment.gotoh.LocalSequenceGotoh;
import huberdp.RDPSolutionTreeOrNode;
import huberdp.Scoring;

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
	public double score(RDPSolutionTreeOrNode node) {
		LocalSequenceGotoh gotoh = new LocalSequenceGotoh(
				-10.0, -2.0,
				bioinfo.alignment.matrices.QuasarMatrix.DAYHOFF_MATRIX);
		
		return ((RDPSolutionTreeOrNode)node.getParent().getParent()).getScore()
				+ gotoh.align(
						node.getProblem().templateSequence,
						node.getProblem().targetSequence
				).getScore();
	}

}

/******************************************************************************
 * "A question that sometimes drives me hazy:                                 *
 *  Am I or are the others crazy?"                                            *
 *     - Albert Einstein (1879 - 1955)                                        *
 ******************************************************************************/
