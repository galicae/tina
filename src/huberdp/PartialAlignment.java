/******************************************************************************
 * huberdp.PartialAlignment.java                                              *
 * Contains the PartialAlignment class which is the definition of an oracle's *
 * partial alignment.                                                                   *
 *                                                                            *
 * This file is best read at line width 80 and tab width 4.                   *
 *                                                                   huberste *
 ******************************************************************************/
package huberdp;

import java.util.LinkedList;

import bioinfo.Sequence;
import bioinfo.alignment.SequenceAlignment;
import bioinfo.proteins.PDBEntry;

/**
 * PartialAlignment is an implementation of the Partial Alignment needed for the
 * RDP data structure.
 * @author huberste
 * @lastchange 2013-02-14
 */
public class PartialAlignment extends RDPProblem {

	/**
	 * start of the partial alignment in the template sequence
	 */
	public int paTemStart;
	/**
	 * end of the partial alignment in the template sequence
	 */
	public int paTemEnd;
	/**
	 * start of the partial alignment in the target sequence
	 */
	public int paTarStart;
	/**
	 * end of the partial alignment in the target sequence
	 */
	public int paTarEnd;
	
	/**
	 * Constructs a new PartialAlignment given the following parameters:
	 * @param templateSequence Sequence of the template
	 * @param templateStructure Structure of the template
	 * @param targetSequence Sequence of the target
	 * @param targetStructure (incomplete) structure of the target (maybe null)
	 * @param partialAlignment partial SequenceAlignment
	 * @param templateStart start of the subproblen (on template)
	 * @param templateEnd end of the subproblem (on template)
	 * @param targetStart start of the subproblem (on target)
	 * @param targetEnd end of the subproblem (on target)
	 * @param paTemStart start of the new partial alignment (template)
	 * @param paTemEnd end of the new partial alignment (template)
	 * @param paTarStart start of the new partial alignment (target)
	 * @param paTarEnd end of the new partial alignment (target)
	 */
	public PartialAlignment(Sequence templateSequence, PDBEntry templateStructure,
			Sequence targetSequence, PDBEntry targetStructure,
			SequenceAlignment partialAlignment,
			int templateStart, int templateEnd, int targetStart, int targetEnd,
			int paTemStart, int paTemEnd, int paTarStart, int paTarEnd) {
		
		super(templateSequence, templateStructure, targetSequence, targetStructure, 
				partialAlignment,
				templateStart, templateEnd, targetStart, targetEnd);
		this.paTemStart = paTemStart;
		this.paTemEnd = paTemEnd;
		this.paTarStart = paTarStart;
		this.paTarEnd = paTarEnd;
	}
	
	/**
	 * Constructs a new PartialAlignment which is the same as the given one. 
	 * @param pa
	 */
	public PartialAlignment(PartialAlignment pa) {
		this(
			pa.templateSequence, pa.templateStructure,
			pa.targetSequence, pa.targetStructure,
			pa.alignment,
			pa.templateStart, pa.templateEnd,
			pa.targetStart, pa.targetEnd,
			pa.paTemStart, pa.paTemEnd,
			pa.paTarStart, pa.paTarEnd
		);
	}

	/**
	 * 
	 * @return the two subproblems that are defined by the partial alignment
	 */
	public LinkedList<RDPProblem> getSubProblems() {
		LinkedList<RDPProblem> results = new LinkedList<RDPProblem>();
		
		// DONE check if this is correct!
		// this seems correct. ~huberste 2013-02-11
		
		int temStart = this.templateStart;
		int temEnd = this.paTemStart-1;
		int tarStart = this.targetStart;
		int tarEnd = this.paTarStart-1;
		
		// DONE check sanity of subproblems! (begin < end, ...)
		if (temStart <= temEnd && tarStart <= tarEnd) {
			results.add (
				new RDPProblem (
					this.templateSequence, this.templateStructure,
					this.targetSequence, this.targetStructure,
					this.alignment,
					temStart, temEnd, tarStart, tarEnd
				)
			);
		}
		
		temStart = this.paTemEnd+1;
		temEnd = this.templateEnd;
		tarStart = this.paTarEnd+1;
		tarEnd = this.targetEnd;
		// DONE check sanity of subproblems! (begin < end, ...)
		if (temStart <= temEnd && tarStart <= tarEnd) {
			results.add (
				new RDPProblem (
					this.templateSequence, this.templateStructure,
					this.targetSequence, this.targetStructure,
					this.alignment,
					temStart, temEnd, tarStart, tarEnd
				)
			);
		}
		
		return results;
	}
	
	/**
	 * toString() function. Mainly for debugging purposes.
	 * @return String representation of the partial alignment.
	 */
	public String toString() {
		String result = alignment.getRowAsString(0)+"\n";
		result += alignment.getRowAsString(1);
		return result;
	}
	
}

/******************************************************************************
 * "A question that sometimes drives me hazy:                                 *
 *  Am I or are the others crazy?"                                            *
 *     - Albert Einstein (1879 - 1955)                                        *
 ******************************************************************************/
