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
		
		LinkedList<PartialAlignment> results = new LinkedList<PartialAlignment>();
		
		// map is the map of the newly generated alignment
		int[][] map = alignment.calcMap();
		// get first aligned character of the template
		
		// calculate paTarStart, paTarEnd, paTemStart, paTemEnd
		// these are the first/last aligned characters of the sequences
		int i = 0;
		// paTemStart is the first aligned character of the template sequence in this subproblem
		int paTemStart = problem.templateStart;
		i=0;
		while (i < map[0].length && map[0][i] < 0) {
			i++;
		}
		if (i < map[0].length) {
			paTemStart += i;
		}
		
		// paTemEnd is the last aligned character of the template sequence in this subproblem
		int paTemEnd = problem.templateEnd;
		i=0;
		while (i < map[0].length && map[0][map[0].length-i-1] < 0) {
			i++;
		}
		if (i < map[0].length) {
			paTemEnd -= i;
		}
		
		// paTarStart is the first aligned character of the target sequence in this subproblem
		int paTarStart = problem.targetStart;
		i=0;
		while (i < map[1].length && map[1][i] < 0) {
			i++;
		}
		if (i < map[1].length) {
			paTarStart += i;
		}
		
		// paTarEnd is the last aligned character of the template sequence in this subproblem
		int paTarEnd = problem.targetEnd;
		i = 0;
		while (i < map[1].length && map[1][map[1].length-i-1] < 0) {
			i++;
		}
		if (i < map[1].length) {
			paTarEnd -= i;
		}
		
		results.add(new PartialAlignment
				(problem.templateSequence, problem.templateStructure,
				problem.targetSequence, problem.targetStructure,
				alignment,
				problem.templateStart, problem.templateEnd,
				problem.targetStart, problem.targetEnd,
				paTemStart, paTemEnd, paTarStart, paTarEnd));
		
//		if (alignment == null) {
//			results = new LinkedList<PartialAlignment>();
//		} else {
//			results = HubeRDP.mergePaA(problem, alignment);
//		}
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