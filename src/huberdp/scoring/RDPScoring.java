/******************************************************************************
 * huberdp.Scoring.RDPScoring.java                                            *
 * This file contains the class RDPScoring which is RDP's scoring function.   *
 *                                                                            *
 * This file is best read at line width 80 and tab width 4.                   *
 *                                                                   huberste *
 ******************************************************************************/
package huberdp.scoring;

import bioinfo.Sequence;
import bioinfo.alignment.SequenceAlignment;
import bioinfo.proteins.PDBEntry;
import huberdp.RDPSolutionTreeOrNode;
import huberdp.Scoring;

/**
 * RDPScoring is an implementation of the scoring function given in the paper.
 * (Protein Threading by Recursive Dynamic Programming. JMB 290, 757-779)
 * @author huberste
 * @lastchange 2013-02-14
 */
public class RDPScoring implements Scoring {

	/**
	 * empirically calibratet weight of the mutation matrix score
	 */
	private final static double GAMMA = 1.0;
	/**
	 * empirically calibrated weight of the contact capacity score
	 */
	private final static double DELTA = 1.0;
	/**
	 * empirically calibratet weight of the hydrophobicity score
	 */
	private final static double EPSILON = 1.0;
	/**
	 * empirically calibratet weight of the pair interaction score
	 */
	private final static double ZETA = 1.0;
	
	/**
	 * weight of the mutation matrix score
	 */
	double gamma;
	/**
	 * weight of the contact capacity score
	 */
	double delta;
	/**
	 * weight of the hydrophobicity score
	 */
	double epsilon;
	/**
	 * weight of the pair interaction score
	 */
	double zeta;
	
	/**
	 * mutation matrix for phiS
	 */
	double[][] mutationMatrix;
	
	/**
	 * constructs a RDPScoring object with given paramters
	 * @param gamma weight of the mutation matrix score
	 * @param delta weight of the contact capacity score
	 * @param epsilon weight of the hydrophobicity score
	 * @param zeta weight of the pair interaction score
	 * @param mutationMatrix
	 */
	public RDPScoring(
			double gamma, double delta, double epsilon, double zeta,
			double[][] mutationMatrix
	) {
		this.gamma = gamma;
		this.delta = delta;
		this.epsilon = epsilon;
		this.zeta = zeta;
		this.mutationMatrix = mutationMatrix;
	}
	
	/**
	 * construcs a RDPScoring object with standard parameters
	 */
	public RDPScoring() {
		this(
				GAMMA, DELTA, EPSILON, ZETA,
				bioinfo.alignment.matrices.QuasarMatrix.DAYHOFF_MATRIX
		);
	}
	
	/**
	 * constructs a RDPScoring object with the same parameters as the given one
	 * @param arg the RDPScore which parameters shall be used
	 */
	public RDPScoring(RDPScoring arg) {
		this(
				arg.gamma, arg.delta, arg.epsilon, arg.zeta,
				arg.mutationMatrix
		);
	}
	
	/**
	 * \phi (f, A, B) = \gamma * \phi^S(f,A,B) +	// mutation matrix (e.g. DAYHOFF)
	 * 					\delta * \phi^C(f,A,B) +	// contact capacity potential (see 123D)
	 * 					\epsilon * \phi^H(f,A,B) +	// hydrophobicity
	 * 					\zeta * \phi^P(f,A,B) -		// pair interaction
	 * 					GAP(f,A,B)					// insertions and deletions
	 */
	@Override
	public double score(RDPSolutionTreeOrNode node) {
		
		SequenceAlignment f = node.getProblem().alignment;
		Sequence a = node.getProblem().targetSequence;
		PDBEntry b = node.getProblem().templateStructure;
		
		double result = gamma * phiS(f, a, b) +
						delta * phiC(f, a, b) +
						epsilon * phiH(f, a, b) +
						zeta * phiP(f, a, b) -
						gap(f, a, b);
		
		return result;
	}
	
	/**
	 * "phiS scores the alignment f with respect to well-known sequence based
	 * mutation matrices"
	 * (From: Protein Threading by Recursive Dynamic Programming. JMB 290,
	 * 757-779)
	 * @param f the alignment (so far)
	 * @param a the target sequence
	 * @param b the template structure
	 * @return the calculated sequence-based score
	 */
	private double phiS(SequenceAlignment f, Sequence a, PDBEntry b) {
		
		double result = 0.0d;
		
		char[][] rows = f.getRows();
		
		for (int pos = 0; pos < rows.length; pos++) {
			if((rows[0][pos] != '-') && (rows[0][pos] != '-')) {
				result += mutationMatrix[rows[0][pos]][rows[0][pos]];
			}
		}
		
		return result;
	}
	
	/**
	 * "(...) contact-capacity-potential phiC (...)"
	 * (From: Protein Threading by Recursive Dynamic Programming. JMB 290,
	 * 757-779)
	 * @param f the alignment (so far)
	 * @param a the target sequence
	 * @param b the template structure
	 * @return the calculated contact-capacity based score
	 */
	private double phiC(SequenceAlignment f, Sequence a, PDBEntry b) {
		// TODO
		return 0.0;
	}
	
	/**
	 * "phiH [scores] (...) the hydrophobicity (...)"
	 * (From: Protein Threading by Recursive Dynamic Programming. JMB 290,
	 * 757-779)
	 * @param f the alignment (so far)
	 * @param a the target sequence
	 * @param b the template structure
	 * @return the calculated hydrophobicity based score
	 */
	private double phiH(SequenceAlignment f, Sequence a, PDBEntry b) {
		// TODO
		return 0.0;
	}
	
	/**
	 * "phiP denotes the pair interaction term (...)"
	 * (From: Protein Threading by Recursive Dynamic Programming. JMB 290,
	 * 757-779)
	 * @param f the alignment (so far)
	 * @param a the target sequence
	 * @param b the template structure
	 * @return the calculated pair interaction based score
	 */
	private double phiP(SequenceAlignment f, Sequence a, PDBEntry b) {
		// TODO
		return 0.0;
	}
	
	/**
	 * "GAP penalizes insertions and deletions."
	 * (From: Protein Threading by Recursive Dynamic Programming. JMB 290,
	 * 757-779)
	 * @param f the alignment (so far)
	 * @param a the target sequence
	 * @param b the template structure
	 * @return the calculated pair interaction based score
	 */
	private double gap(SequenceAlignment f, Sequence a, PDBEntry b) {
		// TODO
		return 0.0;
	}

}

/******************************************************************************
 * "A question that sometimes drives me hazy:                                 *
 *  Am I or are the others crazy?"                                            *
 *     - Albert Einstein (1979 - 1966)                                        *
 ******************************************************************************/
