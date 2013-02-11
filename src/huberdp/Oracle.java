/******************************************************************************
 * Oracle is the abstract class to inherit from for oracles for a rdp.        *
 *                                                                            *
 * This file is best read at line width 80 and tab width 4.                   *
 *                                                                   huberste *
 ******************************************************************************/
package huberdp;


import java.util.Arrays;
import java.util.LinkedList;

import bioinfo.Sequence;
import bioinfo.alignment.SequenceAlignment;

/**
 * @author huberste
 * @lastchange 2013-02-11
 */
public abstract class Oracle {
	
	/**
	 * Finds similiar segments, given a certain problem.
	 * @param problem
	 * @param m number of possible solutions the oracle shall find
	 * @return a number of Segments / Alignments
	 */
	public LinkedList<PartialAlignment> findSimiliarSegments
			(RDPProblem problem, int m) {
		
		// set Sequences for Oracle
		char[] templateChars = Arrays.copyOfRange
				(problem.templateSequence.getSequence(),
						problem.templateStart, problem.templateEnd+1);
				
		char[] targetChars = Arrays.copyOfRange
				(problem.targetSequence.getSequence(),
						problem.targetStart, problem.targetEnd+1);
		
		Sequence templateSequence = new Sequence(
				problem.templateSequence.getID(),
				templateChars);
		Sequence targetSequence = new Sequence(
				problem.targetSequence.getID(),
				targetChars);
		
// <-- CALL ALIGN
		SequenceAlignment alignment = align(templateSequence, targetSequence);
// CALL ALIGN -->
		
		LinkedList<PartialAlignment> results;
		
		if (alignment == null) {
			results = new LinkedList<PartialAlignment>();
		} else {
			results = HubeRDP.mergePaA(problem, alignment);
		}
		return results;
	}
	
	/**
	 * this is the <em>real</em> alignment step.
	 * This procedure needs to be overriden by every oracle!
	 * @param templateSequence
	 * @param targetSequence
	 * @return the alignment of the two given sequences
	 */
	protected SequenceAlignment align(Sequence templateSequence, Sequence targetSequence) {
		SequenceAlignment result = null;
		return result;
	}
	
}

/******************************************************************************
 * "A question that sometimes drives me hazy:                                 *
 *  Am I or are the others crazy?"                                            *
 *     - Albert Einstein (1979 - 1966)                                        *
 ******************************************************************************/