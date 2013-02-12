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
 * @lastchange 2013-02-12
 * 
 */
public class RDPProblem {
	
	/**
	 * Template
	 */
	public Sequence templateSequence;
	public PDBEntry templateStructure;
	
	/**
	 * Target
	 */
	public Sequence targetSequence;
	public PDBEntry targetStructure;
	
	/**
	 * SequenceAlignment so far
	 */
	public SequenceAlignment alignment;
	
	/**
	 * start and end of the given (sub-)problem.
	 */
	public int templateStart;
	public int templateEnd;
	public int targetStart;
	public int targetEnd;
	
	/**
	 * Constructs an RDPProblem
	 * @param templateSequence
	 * @param templateStructure
	 * @param targetSequence
	 * @param targetStructure
	 * @param partialAlignment
	 * @param templateStart Start of the (sub-) problem on template side
	 * @param templateEnd End of the (sub-) problem on template side
	 * @param targetStart Start of the (sub-) problem on target side
	 * @param targetEnd End of the (sub-) problem on target side
	 */
	public RDPProblem(
			Sequence templateSequence, PDBEntry templateStructure,
			Sequence targetSequence, PDBEntry targetStructure,
			SequenceAlignment partialAlignment,
			int templateStart, int templateEnd,
			int targetStart, int targetEnd) {
		this.templateSequence = templateSequence;
		this.templateStructure = templateStructure;
		this.targetSequence = targetSequence;
		this.targetStructure = targetStructure;
		this.alignment = partialAlignment;
		this.templateStart = templateStart;
		this.templateEnd = templateEnd;
		this.targetStart = targetStart;
		this.targetEnd = targetEnd;
	}
	
	/**
	 * Constructs a new RDPProblem which is a copy of the given one.
	 * @param problem
	 */
	public RDPProblem(RDPProblem problem) {
		this.templateSequence = problem.templateSequence;
		this.templateStructure = problem.templateStructure;
		this.targetSequence = problem.targetSequence;
		this.targetStructure = problem.targetStructure;
		this.alignment = problem.alignment;
		this.templateStart = problem.templateStart;
		this.templateEnd = problem.templateEnd;
		this.targetStart = problem.targetStart;
		this.targetEnd = problem.targetEnd;
	}
	
	/**
	 * toString() function. mainly for debugging.
	 */
	public String toString() {
		String result = "";
		result += templateSequence.getSequenceAsString().substring(templateStart, templateEnd+1)+ "\n";
		result += targetSequence.getSequenceAsString().substring(targetStart, targetEnd+1);
		return result;
	}
}

/******************************************************************************
 * "A question that sometimes drives me hazy:                                 *
 *  Am I or are the others crazy?"                                            *
 *     - Albert Einstein (1979 - 1966)                                        *
 ******************************************************************************/