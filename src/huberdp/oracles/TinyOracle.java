/******************************************************************************
 * huberdp.oracles.TinyOracle.java                                            *
 * Cotains the class TinyOracle which is an oracle for the RDP featuring a    *
 * local Gotoh.                                                               *
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
 * TinyOracle is an Oracle for HubeRDP that features an local Gotoh as sequence
 * alignment prediction.
 * @author huberste
 * @lastchange 2013-02-14
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

}

/******************************************************************************
 * "A question that sometimes drives me hazy:                                 *
 *  Am I or are the others crazy?"                                            *
 *     - Albert Einstein (1979 - 1966)                                        *
 ******************************************************************************/
