/******************************************************************************
 * RDPProblem is the problem definition for RDP.                              *
 *                                                                            *
 * This file is best read at line width 80 and tab width 4.                   *
 *                                                                   huberste *
 ******************************************************************************/
package huberdp;

import bioinfo.Sequence;
import bioinfo.alignment.SequenceAlignment;
import bioinfo.proteins.PDBEntry;

/**
 * @author huberste
 * @lastchange 2013-01-26
 * 
 */
public class RDPProblem {
	
	/**
	 * Target
	 */
	public Sequence targetSequence;
	public PDBEntry targetStructure;
	
	/**
	 * Template
	 */
	public Sequence templateSequence;
	public PDBEntry templateStructure;
	
	/**
	 * SequenceAlignment so far
	 */
	public SequenceAlignment alignment;
	
	/**
	 * start of the given (sub-)problem.
	 */
	public int targetStart;
	public int templateStart;
	/**
	 * end of the given (sub-)problem.
	 */
	public int targetEnd;
	public int templateEnd;
	
	/**
	 * Constructs an RDPProblem
	 */
	public RDPProblem(
			Sequence targetSequence, PDBEntry targetStructure,
			Sequence templateSequence, PDBEntry templateStructure,
			SequenceAlignment partialAlignment,
			int targetStart, int targetEnd, int templateStart, int templateEnd) {
		this.targetSequence = targetSequence;
		this.targetStructure = targetStructure;
		this.templateSequence = templateSequence;
		this.templateStructure = templateStructure;
		this.alignment = partialAlignment;
		this.targetStart = targetStart;
		this.targetEnd = targetEnd;
		this.templateStart = templateStart;
		this.templateEnd = templateEnd;
		
	}
}

/******************************************************************************
 * "A question that sometimes drives me hazy:                                 *
 *  Am I or are the others crazy?"                                            *
 *     - Albert Einstein (1979 - 1966)                                        *
 ******************************************************************************/