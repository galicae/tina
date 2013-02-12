/******************************************************************************
 * TreeAlignment is part of the problem solution of RDP.                      *
 *                                                                            *
 * This file is best read at line width 80 and tab width 4.                   *
 *                                                                   huberste *
 ******************************************************************************/
package huberdp;

import bioinfo.Sequence;
import bioinfo.alignment.Alignment;
import bioinfo.proteins.PDBEntry;

/**
 * This class is the solution of an RDPProblem.
 * @author huberste
 * @lastchange 2013-02-11
 */
public class TreeAlignment{
	
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
	 * the tree alignment
	 */
	public Alignment alignment;
	
	/**
	 * constructs a  new RDPSolution
	 * @param targetSequence
	 * @param targetStructure
	 */
	public TreeAlignment(Sequence templateSequence, PDBEntry templateStructure,
			Sequence targetSequence, PDBEntry targetStructure,
			Alignment alignment) {
		this.templateSequence = templateSequence;
		this.templateStructure = templateStructure;
		this.targetSequence = targetSequence;
		this.targetStructure = targetStructure;
		this.alignment = alignment;
	}
	
	public String toString() {
		String result = "RDPSolution\n";
		result += alignment.toStringVerbose();
		return result;
	}
	
}

/******************************************************************************
 * "A question that sometimes drives me hazy:                                 *
 *  Am I or are the others crazy?"                                            *
 *     - Albert Einstein (1979 - 1966)                                        *
 ******************************************************************************/