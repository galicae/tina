/******************************************************************************
 * huberdp.TreeAlignment.java                                                 *
 *                                                                            *
 * Contains the class TreeAlignment which is part of the problem solution of  *
 * RDP.                                                                       *
 *                                                                            *
 * This file is best read at line width 80 and tab width 4.                   *
 *                                                                   huberste *
 ******************************************************************************/
package huberdp;

import bioinfo.alignment.Threading;

/**
 * TreeAlignment represents the TreeAlignment structure of RDP.
 * 
 * @author huberste
 * @lastchange 2013-02-19
 */
public class TreeAlignment {

	private Threading threading;

	/**
	 * Standard Constructor for a new TreeAlignment
	 * 
	 * @param threading
	 */
	public TreeAlignment(Threading threading) {
		this.threading = threading;
	}

	public Threading getThreading() {
		return threading;
	}

	/**
	 * toString() method. mostly for debugging.
	 */
	public String toString() {
		String result = "TA\n";
		result += threading.toStringVerbose();
		return result;
	}

}

/******************************************************************************
 * "A question that sometimes drives me hazy:                                 *
 *  Am I or are the others crazy?"                                            *
 * - Albert Einstein (1879 - 1955)                                            *
 ******************************************************************************/
