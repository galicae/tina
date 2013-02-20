/******************************************************************************
 * huberdp.RDPProblem.java                                                    *
 * Contains the class RDPProblem which is the (sub-)problem definition for an *
 * OR node of RDP.                                                            *
 *                                                                            *
 * This file is best read at line width 80 and tab width 4.                   *
 *                                                                   huberste *
 ******************************************************************************/
package huberdp;

import bioinfo.alignment.Threading;

/**
 * RDPProblem is an implementation of an (sub-)problem defined by an OR node.
 * 
 * @author huberste
 * @lastchange 2013-02-14
 */
public class RDPProblem {

	/**
	 * current Threading
	 */
	private Threading threading;
	/**
	 * Start of the currend (sub-)problem, inclusive
	 */
	private int problemStart;
	/**
	 * End of the current (sub-)problem, inclusive
	 */
	private int problemEnd;

	/**
	 * Constructs an RDPProblem
	 * 
	 * @param threading
	 *            the Threading that defines this (sub-)problem
	 * @param problemStart
	 *            start of the (sub-)problem
	 * @param problemEnd
	 *            end of the (sub-)problem
	 */
	public RDPProblem(Threading threading, int problemStart, int problemEnd) {
		this.threading = threading;
		this.problemStart = problemStart;
		this.problemEnd = problemEnd;
	}

	/**
	 * 
	 * @return the Threading
	 */
	public Threading getThreading() {
		return threading;
	}

	/**
	 * 
	 * @return the problemStart value
	 */
	public int getProblemStart() {
		return problemStart;
	}

	/**
	 * 
	 * @return the problemEnd value
	 */
	public int getProblemEnd() {
		return problemEnd;
	}

	/**
	 * toString() function. mainly for debugging.
	 */
	public String toString() {
		String[] result = new String[2];
		String[] rows = threading.getRowsAsString();
		result[0] = threading.getComponent(0).getID() + ": ";
		result[1] = threading.getComponent(1).getID() + ": ";
		for (int i = 0; i < threading.length(); i++) {
			if (problemStart == i) {
				result[0] += ">";
				result[1] += ">";
			}
			result[0] += rows[0].charAt(i);
			result[1] += rows[1].charAt(i);
			if (problemEnd == i) {
				result[0] += "<";
				result[1] += "<";
			}

		}
		return result[0] + "\n" + result[1] + "\n";
	}
}

/******************************************************************************
 * "A question that sometimes drives me hazy: * Am I or are the others crazy?" *
 * - Albert Einstein (1879 - 1955) *
 ******************************************************************************/
