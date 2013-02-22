/******************************************************************************
 * huberdp.Scoring.java                                                       *
 *                                                                            *
 * Contains the interface Scoring which is to be implemented by every Scoring *
 * class for scoring Threadings                                               *
 *                                                                            *
 * This file is best read at line width 80 and tab width 4.                   *
 *                                                                   huberste *
 ******************************************************************************/
package huberdp;

import bioinfo.alignment.Threading;

/**
 * Scoring is an interface for scoring Threadings
 * 
 * Possible Scoring features: * unclosable gaps * large gaps * Hydrophobicity *
 * ContactCapacityPotentials * Pair Contact Potentials * Sequence Scoring
 * Matrices
 * 
 * @author huberste
 * @lastchange 2013-02-22
 */
public interface Scoring {

	/**
	 * scores a complete Threading
	 * 
	 * @param threading
	 * @return
	 */
	public double score(Threading threading);

	/**
	 * scores a single match
	 * 
	 * @param threading
	 * @param m
	 *            position in template
	 * @param n
	 *            position in target
	 * @return
	 */
	public double getScore(Threading t, int m, int n);

	/**
	 * scores a single insertion. Score should be nagative.
	 * 
	 * @param t
	 *            threading
	 * @param n
	 *            position of the inserted AminoAcid in the target
	 * @return
	 */
	public double getInsertionScore(Threading t, int n);

	/**
	 * scores a single deletion. Score should be negative.
	 * 
	 * @param t
	 *            threading
	 * @param m
	 *            position of the deleted AminoAcid in the template
	 * @return
	 */
	public double getDeletionScore(Threading t, int m);

}

/******************************************************************************
 * "A question that sometimes drives me hazy:                                 *
 *  Am I or are the others crazy?"                                            *
 * - Albert Einstein (1879 - 1955)                                            *
 ******************************************************************************/
