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
		
		if (alignment == null) {
			// TODO: cannot be further aligned by this oracle.
		} else {
		
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
		
		// calculate new partial alignment (pa)
		SequenceAlignment pa = null;
		if (problem.alignment == null) {
			// no pa was calculated yet, this means that this is the first call.
			pa = alignment;
		} else {
			// TODO check this code!
			// merge problem.alignment with alignment
			map = problem.alignment.calcMap();
			char[] temRow = problem.alignment.getRow(0); // templateRow
			char[] tarRow = problem.alignment.getRow(1); // targetRow
			
			int paTemStartInd = 0; // position in temRow, later first aligned character
			int pos = 0; // position in problem.templateSequence
			while (pos < problem.templateStart) {
				if (temRow[paTemStartInd] != '-') {
					pos++;
				}
				paTemStartInd++;
			}
			
			int paTarStartInd = 0; // position in tarRow, later first aligned character
			pos = 0; // position in problem.targetSequence
			while (pos < problem.targetStart) {
				if (tarRow[paTarStartInd] != '-') {
					pos++;
				}
				paTarStartInd++;
			}
			
			int paTemEndInd = temRow.length - 1; // position in temRow, later last aligned character
			pos = problem.templateSequence.length() - 1; // position in problem.templateSequence
			while (pos > problem.templateEnd) {
				if (temRow[paTemEndInd] != '-') {
					pos--;
				}
				paTemEndInd--;
			}
			
			int paTarEndInd = tarRow.length - 1; // position in tarRow, later last aligned character
			pos = problem.targetSequence.length() - 1;  // position in problem.targetSequence
			while (pos > problem.targetEnd) {
				if (tarRow[paTarEndInd] != '-') {
					pos--;
				}
				paTarEndInd--;
			}
			
			// Error handling 143
			if (paTemStartInd != paTarStartInd) {
				System.err.println("Error 143 in Oracle: Alignment lengths don't match.");
			}
			// Error handling 147
			if (paTemEndInd != paTarEndInd) {
				System.err.println("Error 147 in Oracle: Alignment lengths don't match.");
			}
			
			char[][] newRows = new char[2][];
			// paTemStartInd = length of old alignment before new alignment
			// alignment.length = length of new alignment
			// (problem.alignment.length() - (paTemEndInd+1)) = length of old alignment behind new alignment
			newRows[0] = new char[paTemStartInd + alignment.length() + (problem.alignment.length() - (paTemEndInd+1))];
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
			// copy end of old alignment
			for (pos = 0; pos < problem.alignment.length() - (paTemEndInd + 1); pos++) {
				newRows[0][paTemStartInd + alignment.length() + pos]=temRow[(paTemEndInd + 1) + pos];
			}
			for (pos = 0; pos < problem.alignment.length() - (paTarEndInd + 1); pos++) {
				newRows[1][paTarStartInd + alignment.length() + pos]=tarRow[(paTarEndInd + 1) + pos];
			}
			
			// merge scores. Here simply adding the scores will do fine.
			double newScore = problem.alignment.getScore()+alignment.getScore();
			
			pa = new SequenceAlignment(
					problem.templateSequence, problem.targetSequence, newRows, newScore);
			
		}
		
		results.add(new PartialAlignment
				(problem.templateSequence, problem.templateStructure,
				problem.targetSequence, problem.targetStructure,
				pa,
				problem.templateStart, problem.templateEnd,
				problem.targetStart, problem.targetEnd,
				paTemStart, paTemEnd, paTarStart, paTarEnd));
		
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