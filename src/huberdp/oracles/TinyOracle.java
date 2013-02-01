/******************************************************************************
 * TinyOracle is an oracle for the RDP featuring an local Gotoh.              *
 *                                                                            *
 * This file is best read at line width 80 and tab width 4.                   *
 *                                                                   huberste *
 ******************************************************************************/
package huberdp.oracles;

import java.util.Arrays;
import java.util.LinkedList;

import bioinfo.Sequence;
import bioinfo.alignment.SequenceAlignment;
import bioinfo.alignment.gotoh.LocalSequenceGotoh;

import huberdp.Oracle;
import huberdp.PartialAlignment;
import huberdp.RDPProblem;

/**
 * @author huberste
 * @lastchange 2013-01-29
 *
 */
public class TinyOracle implements Oracle {

	/* (non-Javadoc)
	 * @see huberdp.Oracle#findSimiliarSegments(huberdp.RDPProblem, int)
	 */
	@Override
	public LinkedList<PartialAlignment> findSimiliarSegments
			(RDPProblem problem, int m) {
		
		// set Sequences for Gotoh
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
		
		// Create new Gotoh
		LocalSequenceGotoh gotoh = new LocalSequenceGotoh(
				-10.0, -2.0,
				bioinfo.alignment.matrices.QuasarMatrix.DAYHOFF_MATRIX);
		SequenceAlignment alignment = gotoh.align
				(templateSequence, targetSequence);
		
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
		while (i >= 0 && map[0][map[0].length-i-1] < 0) {
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
			// TODO merge problem.alignment with alignment
			map = problem.alignment.calcMap();
			problem.alignment.getRow(0);
			problem.alignment.getRow(1);
		}
		
		
		char[][] newRows = new char[2][];
		
		// merge scores. Here simply adding the scores will be fine.
		double newScore = problem.alignment.getScore()+alignment.getScore();
		
		pa = new SequenceAlignment(
				problem.templateSequence, problem.targetSequence, newRows, newScore);
		

		
		results.add(new PartialAlignment
				(problem.templateSequence, problem.targetStructure,
				problem.templateSequence, problem.templateStructure,
				pa,
				problem.targetStart, problem.targetEnd,
				problem.templateStart, problem.templateEnd,
				paTarStart, paTarEnd, paTemStart, paTemEnd));
		
		return results;
	}

}
