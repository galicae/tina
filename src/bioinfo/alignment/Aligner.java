/******************************************************************************
 * bioinfo.alignment.Aligner                                                  *
 *                                                                            *
 * This file is best read at line width 80 and tab width 4.                   *
 ******************************************************************************/
package bioinfo.alignment;

/**
 * @author gobi_4
 * @date November 25, 2012
 */
public interface Aligner {
	
	
	/**
	 * This method should align the two given sequences.
	 * @param sequence1
	 * @param sequence2
	 * @return the Alignment of the two given sequences
	 */
	public Alignment align(Alignable sequence1, Alignable sequence2);
	
	/**
	 * This method should check if the score in the given Alignment is correct.
	 * @param arg
	 * @return true if score is correct, false else.
	 */
	public abstract boolean check(Alignment arg);
	
}