/**
 * 
 */
package huberdp.scoring;

import bioinfo.Sequence;
import bioinfo.alignment.Alignment;
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
	
	public RDPScoring() {
		this.gamma = GAMMA;
		this.delta = DELTA;
		this.epsilon = EPSILON;
		this.zeta = ZETA;
	}
	
	public RDPScoring(double gamma, double delta, double epsilon, double zeta) {
		this.gamma = gamma;
		this.delta = delta;
		this.epsilon = epsilon;
		this.zeta = zeta;
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
		
		Alignment f = node.getProblem().alignment;
		Sequence a = node.getProblem().targetSequence;
		PDBEntry b = node.getProblem().templateStructure;
		
		double result = gamma * phiS(f, a, b) +
						delta * phiC(f, a, b) +
						epsilon * phiH(f, a, b) +
						zeta * phiP(f, a, b) -
						gap(f, a, b);
		
		return result;
	}
	
	private double phiS(Alignment f, Sequence a, PDBEntry b) {
		// TODO
		return 0.0;
	}
	
	private double phiC(Alignment f, Sequence a, PDBEntry b) {
		// TODO
		return 0.0;
	}
	
	private double phiH(Alignment f, Sequence a, PDBEntry b) {
		// TODO
		return 0.0;
	}
	
	private double phiP(Alignment f, Sequence a, PDBEntry b) {
		// TODO
		return 0.0;
	}
	
	private double gap(Alignment f, Sequence a, PDBEntry b) {
		// TODO
		return 0.0;
	}

}
