/******************************************************************************
 * huberdp.oracles.TinyOracle.java                                            *
 * Cotains the class TinyOracle which is an oracle for the RDP featuring a    *
 * local Gotoh.                                                               *
 *                                                                            *
 * This file is best read at line width 80 and tab width 4.                   *
 *                                                                   huberste *
 ******************************************************************************/
package huberdp.oracles;

import java.util.LinkedList;

import bioinfo.Sequence;
import bioinfo.alignment.SequenceAlignment;
import bioinfo.alignment.gotoh.LocalSequenceGotoh;

import huberdp.Oracle;
import huberdp.PartialAlignment;
import huberdp.RDPProblem;

/**
 * TinyOracle is an Oracle for HubeRDP that features an local Gotoh as sequence
 * alignment prediction.
 * @author huberste
 * @lastchange 2013-02-14
 */
public class TinyOracle implements Oracle {

	/**
	 * this is the _real_ alignment step.
	 * @param templateSequence
	 * @param targetSequence
	 * @return the alignment of the two given sequences
	 */
	protected SequenceAlignment align(Sequence templateSequence, Sequence targetSequence) {
				
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

	@Override
	public LinkedList<PartialAlignment> findSimiliarSegments(
			RDPProblem problem, int m) {
		
		// allocate result
		LinkedList<PartialAlignment> results = new LinkedList<PartialAlignment>();

		// set Sequences for Oracle
		String[] rows = problem.getThreading().getRowsAsString();
		String template = "";
		String target = "";
		for (int i = problem.getProblemStart(); i < problem.getProblemEnd(); i++) {
			if (rows[0].charAt(i) != '-') {
				template += rows[0].charAt(i);
			}
			if (rows[1].charAt(i) != '-') {
				target += rows[1].charAt(i);
			}
		}
		
		Sequence templateSequence = new Sequence(problem.getThreading().getStructure().getID(), template);
		Sequence   targetSequence = new Sequence(problem.getThreading().getSequence().getId(), target);
		results.add(new PartialAlignment(problem, align(templateSequence, targetSequence)));
		return results;
	}

}

/******************************************************************************
 * "A question that sometimes drives me hazy:                                 *
 *  Am I or are the others crazy?"                                            *
 *     - Albert Einstein (1879 - 1955)                                        *
 ******************************************************************************/
