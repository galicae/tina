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
						problem.templateStart, problem.templateEnd);
				
		char[] targetChars = Arrays.copyOfRange
				(problem.targetSequence.getSequence(),
						problem.targetStart, problem.targetEnd);
		
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
		
		if (alignment == null) {
			// TODO: cannot be further aligned by this oracle.
		} else {
		
		int[][] map = alignment.calcMap();
		// get first aligned character of the template
		
		// calculate paTarStart, paTarEnd, paTemStart, paTemEnd
		int i = 0;
		int paTemStart = problem.templateStart;
		i=0;
		while (i < map[0].length && map[0][i] < 0) {
			i++;
		}
		if (i < map[0].length) {
			paTemStart += i;
		}
		
		int paTemEnd = problem.templateEnd;
		i=0;
		while (i < map[0].length && map[0][map[0].length-i-1] < 0) {
			i++;
		}
		if (i < map[0].length) {
			paTemEnd -= i;
		}
		
		int paTarStart = problem.targetStart;
		i=0;
		while (i < map[1].length && map[1][i] < 0) {
			i++;
		}
		if (i < map[1].length) {
			paTarStart += i;
		}
		
		int paTarEnd = problem.targetEnd;
		i = 0;
		while (i < map[1].length && map[1][map[1].length-i-1] < 0) {
			i++;
		}
		if (i < map[1].length) {
			paTarEnd -= i;
		}
		
		// calculate new partial alignment
		SequenceAlignment pa = null;
		if (problem.alignment == null) {
			pa = alignment;
		} else {
			// merge problem.alignment with alignment
			map = problem.alignment.calcMap();
			char[] temRow = problem.alignment.getRow(0);
			char[] tarRow = problem.alignment.getRow(1);
			int paTemStartInd = 0;
			int pos = 0;
			while (pos < paTemStart) {
				if (temRow[paTemStartInd] != '-') {
					pos++;
				}
				paTemStartInd ++;
			}
			int paTarStartInd = 0;
			pos = 0;
			while (pos < paTarStart) {
				if (tarRow[paTarStartInd] != '-') {
					pos++;
				}
				paTarStartInd ++;
			}
			int paTemEndInd = 0;
			pos = 0;
			while (pos < paTemEnd) {
				if (temRow[paTemEndInd] != '-') {
					pos++;
				}
				paTemEndInd ++;
			}
			int paTarEndInd = 0;
			pos = 0;
			while (pos < paTarEnd) {
				if (tarRow[paTarEndInd] != '-') {
					pos++;
				}
				paTarEndInd ++;
			}
			// Error handling 141
			if (paTemStartInd != paTarStartInd) {
				System.err.println("Error 141 in TinyOracle: Alignment lengths don't match.");
			}
			// Error handling 145
			if (paTemEndInd != paTarEndInd) {
				System.err.println("Error 145 in TinyOracle: Alignment lengths don't match.");
			}
			
			char[][] newRows = new char[2][];
			newRows[0] = new char[(paTemStartInd - 1) + alignment.length() + (problem.alignment.length() - paTemEndInd)];
			newRows[1] = new char[newRows[0].length];
			// make new Alignment:
			// copy beginning of old alignment
			for (pos = 0; pos < paTemStartInd; pos++) {
				newRows[0][pos]=temRow[pos];
			}
			for (pos = 0; pos < paTarStartInd; pos++) {
				newRows[1][pos]=tarRow[pos];
			}
			// copy new aligned
			for (pos = 0; pos < alignment.length(); pos++) {
				newRows[0][paTemStartInd + pos]=alignment.getRow(0)[pos];
			}
			for (pos = 0; pos < alignment.length(); pos++) {
				newRows[1][paTarStartInd + pos]=alignment.getRow(1)[pos];
			}
			// copy new aligned
			for (pos = 0; pos < problem.alignment.length() - paTemEndInd; pos++) {
				newRows[0][paTemStartInd + alignment.length() + pos]=temRow[paTemEndInd + pos];
			}
			for (pos = 0; pos < problem.alignment.length() - paTarEndInd; pos++) {
				newRows[1][paTarStartInd + alignment.length() + pos]=tarRow[paTarEndInd + pos];
			}
			
			// merge scores. Here simply adding the scores will be fine.
			double newScore = problem.alignment.getScore()+alignment.getScore();
			
			pa = new SequenceAlignment(
					problem.templateSequence, problem.targetSequence, newRows, newScore);
			
		}
		
		results.add(new PartialAlignment
				(problem.templateSequence, problem.targetStructure,
				problem.templateSequence, problem.templateStructure,
				pa,
				problem.targetStart, problem.targetEnd,
				problem.templateStart, problem.templateEnd,
				paTarStart, paTarEnd, paTemStart, paTemEnd));
		
		}
		
		return results;
	}
	
	/**
	 * this is the <em>real</em> alignment step.
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