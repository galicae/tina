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
		// TODO
		
		// set Sequences for Gotoh
		char[] targetChars = Arrays.copyOfRange
				(problem.targetSequence.getSequence(),
						problem.targetStart, problem.targetEnd+1);
		
		
		char[] templateChars = Arrays.copyOfRange
				(problem.templateSequence.getSequence(),
						problem.templateStart, problem.templateEnd+1);
		
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
		// TODO
		
		LinkedList<RDPProblem> results = new LinkedList<RDPProblem>();
		results.add(null);
		return results;
	}

}
