/******************************************************************************
 * huberdp.scoring.SimpleScoring.java                                         *
 *                                                                            *
 * Contains the class SimpleScoring that provides a simple scoring function   *
 * for OR nodes.                                                              *
 *                                                                            *
 * This file is best read at line width 80 and tab width 4.                   *
 *                                                                   huberste *
 ******************************************************************************/
package huberdp.scoring;

import huberdp.Scoring;
import bioinfo.alignment.Threading;
import bioinfo.alignment.matrices.QuasarMatrix;

/**
 * SimpleScoring provides a simple scoring function for OR nodes.
 * 
 * @author huberste
 * @lastchange 2013-02-14
 */
public class SimpleScoring implements Scoring {

	private double[][] dayhoff = QuasarMatrix.DAYHOFF_MATRIX;

	private boolean gap = false;
	double lastgap = 0.0;

	/*
	 * (non-Javadoc)
	 * 
	 * @see huberdp.Scoring#score(huberdp.RDPSolutionTreeOrNode)
	 */
	@Override
	public double score(Threading threading) {
		double result = 0;
		char[][] rows = threading.getRowsAsCharArray();
		boolean first = false;
		boolean gap = false;
		double lastgap = 0.0;
		for (int i = 0; i < threading.length(); i++) {
			if (rows[0][i] == '-') { // insertion
				if (first) { // first Gap is through
					if (gap) {
						lastgap -= 2.0;
					} else {
						gap = true;
						lastgap = 10;
					}
				}
			} else if (rows[1][i] == '-') { // deletion
				if (first) { // first Gap is through
					if (gap) { // gap already opened
						lastgap -= 2.0;
					} else {
						gap = true;
						lastgap = 10;
					}
				}
			} else { // match
				first = true;
				gap = false;
				result += lastgap;
				result += dayhoff[rows[0][i] - 65][rows[1][i] - 65];
			}
		}
		return result;
	}

	@Override
	public double getScore(Threading t, int m, int n) {
		gap = false;
		double result = lastgap;
		result += dayhoff[t.getStructure().getAminoAcid(m).getName()
				.getNumber()][t.getSequence().getComp(n) - 65];
		lastgap = 0;
		return result;
	}

	@Override
	public double getInsertionScore(Threading t, int n) {
		return getGapScore();
	}

	@Override
	public double getDeletionScore(Threading t, int m) {
		return getGapScore();
	}
	
	private double getGapScore() {
		// check if there comes a match after this gap isn't needed because
		// local
		if (gap) {
			return -2.0;
		} else {
			gap = true;
			return -10;
		}
	}

}

/******************************************************************************
 * "A question that sometimes drives me hazy:                                 *
 *  Am I or are the others crazy?"                                            *
 * - Albert Einstein (1879 - 1955)                                            *
 ******************************************************************************/
