/******************************************************************************
 * TinyOracle is an oracle for the RDP featuring an local Gotoh.              *
 *                                                                            *
 * This file is best read at line width 80 and tab width 4.                   *
 *                                                                   huberste *
 ******************************************************************************/
package huberdp.oracles;

import bioinfo.Sequence;
import bioinfo.alignment.SequenceAlignment;
import bioinfo.alignment.gotoh.LocalSequenceGotoh;

import huberdp.Oracle;

/**
 * @author huberste
 * @lastchange 2013-02-11
 *
 */
public class TinyOracle extends Oracle {

	/**
	 * this is the _real_ alignment step.
	 * @param templateSequence
	 * @param targetSequence
	 * @return the alignment of the two given sequences
	 */
	@Override
	protected SequenceAlignment align(Sequence templateSequence, Sequence targetSequence) {
		// TODO check if correct!
				
		// Create new Gotoh
		LocalSequenceGotoh gotoh = new LocalSequenceGotoh(
				-10.0, -2.0,
				bioinfo.alignment.matrices.QuasarMatrix.DAYHOFF_MATRIX);
		// start debugging
//				System.out.println(templateSequence + "\n" + targetSequence);
		// end debugging
		SequenceAlignment result = gotoh.align
				(templateSequence, targetSequence);
		// begin debugging
//		System.out.println("Gotoh Alignment worked fine:\n"+result.toStringVerbose());
		// end debugging
		
		return result;
	}
	
	/*
	public LinkedList<PartialAlignment> findSimiliarSegments
			(RDPProblem problem, int m) {
		
		
		LinkedList<PartialAlignment> results = new LinkedList<PartialAlignment>();
		
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
		
		return results;
	}
	*/

}
