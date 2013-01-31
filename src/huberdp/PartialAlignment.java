/******************************************************************************
 * PartialAlignment is the definition of an oracle's segment.                 *
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
 * @lastchange 2012-01-29
 *
 */
public class PartialAlignment extends RDPProblem {

	public int paTarStart;
	public int paTarEnd;
	public int paTemStart;
	public int paTemEnd;
	
	/**
	 * Constructs a new PartialAlignment given the following parameters:
	 * @param targetSequence Sequence of the target
	 * @param targetStructure (incomplete) structure of the target (maybe null)
	 * @param templateSequence Sequence of the template
	 * @param templateStructure Structure of the template
	 * @param partialAlignment partial SequenceAlignment
	 * @param targetStart start of the subproblem (on target)
	 * @param targetEnd end of the subproblem (on target)
	 * @param templateStart start of the subproblen (on template)
	 * @param templateEnd end of the subproblem (on template)
	 * @param paTarStart start of the new partial alignment (target)
	 * @param paTarEnd end of the new partial alignment (target)
	 * @param paTemStart start of the new partial alignment (template)
	 * @param paTemEnd end of the new partial alignment (template)
	 */
	public PartialAlignment(Sequence targetSequence, PDBEntry targetStructure,
			Sequence templateSequence, PDBEntry templateStructure,
			SequenceAlignment partialAlignment,
			int targetStart, int targetEnd, int templateStart, int templateEnd,
			int paTarStart, int paTarEnd, int paTemStart, int paTemEnd) {
		super(targetSequence, targetStructure, templateSequence, templateStructure,
				partialAlignment, targetStart, targetEnd, templateStart, templateEnd);
		this.paTarStart = paTarStart;
		this.paTarEnd = paTarEnd;
		this.paTemStart = paTemStart;
		this.paTemEnd = paTemEnd;
	}

	/**
	 * 
	 * @return the two subproblems that are defined by the partial alignment
	 */
	public RDPProblem[] getSubProblems() {
		RDPProblem[] result = new RDPProblem[2];
		
		// TODO check if this is correct!
		int tarStart = this.targetStart;
		int tarEnd = this.paTarStart-1;
		int temStart = this.templateStart;
		int temEnd = this.paTemStart-1;
		result[0] = new RDPProblem
				(this.targetSequence, this.targetStructure,
						this.templateSequence, this.templateStructure,
						this.alignment,
						tarStart, tarEnd, temStart, temEnd);
		
		tarStart = this.paTarEnd+1;
		tarEnd = this.targetEnd;
		temStart = this.paTemEnd+1;
		temEnd = this.templateEnd;
		result[1] = new RDPProblem
				(this.targetSequence, this.targetStructure,
						this.templateSequence, this.templateStructure,
						this.alignment,
						tarStart, tarEnd, temStart, temEnd);
		
		return result;		
	}
	
}
