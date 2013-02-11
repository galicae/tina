/******************************************************************************
 * TreeAlignment is part of the problem solution of RDP.                      *
 *                                                                            *
 * This file is best read at line width 80 and tab width 4.                   *
 *                                                                   huberste *
 ******************************************************************************/
package huberdp;

import bioinfo.Sequence;
import bioinfo.alignment.SequenceAlignment;
import bioinfo.proteins.PDBEntry;

/**
 * I have no Idea what I'm doing.
 * @author huberste
 * @lastchange 2013-02-11
 */
public class RDPSolution{
	
	/**
	 * Target == solution
	 */
	public Sequence targetSequence;
	public PDBEntry targetStructure;
	
	/**
	 * 
	 * @param templateSequence
	 * @param templateStructure
	 * @param targetSequence
	 * @param targetStructure
	 * @param solution
	 */
	public RDPSolution(Sequence templateSequence, PDBEntry templateStructure,
			Sequence targetSequence, PDBEntry targetStructure,
			SequenceAlignment solution) {
		this.targetSequence = targetSequence;
		this.targetStructure = targetStructure;
	}

	
	
}

/******************************************************************************
 * "A question that sometimes drives me hazy:                                 *
 *  Am I or are the others crazy?"                                            *
 *     - Albert Einstein (1979 - 1966)                                        *
 ******************************************************************************/