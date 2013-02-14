/******************************************************************************
 * huberdp.Scoring.RDPScoring.java                                            *
 * This file contains the class RDPScoring which is RDP's scoring function.   *
 *                                                                            *
 * This file is best read at line width 80 and tab width 4.                   *
 *                                                                   huberste *
 ******************************************************************************/
package huberdp.scoring;

import java.util.HashMap;
import java.util.Set;

import bioinfo.Sequence;
import bioinfo.alignment.SequenceAlignment;
import bioinfo.energy.potential.voronoi.VoroPPWrap;
import bioinfo.energy.potential.voronoi.VoroPrepType;
import bioinfo.energy.potential.voronoi.VoronoiData;
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
	 * vorobin absolute path to location of voro++ binary (dont have one? look at ./tools/voro++)
	 */
	String vorobin;
	
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
			double[][] mutationMatrix, String vorobin
	) {
		this.gamma = gamma;
		this.delta = delta;
		this.epsilon = epsilon;
		this.zeta = zeta;
		this.mutationMatrix = mutationMatrix;
		this.vorobin = vorobin;
	}
	
	/**
	 * construcs a RDPScoring object with standard parameters
	 */
	public RDPScoring() {
		this(
				GAMMA, DELTA, EPSILON, ZETA,
				bioinfo.alignment.matrices.QuasarMatrix.DAYHOFF_MATRIX,null
		);
	}
	
	/**
	 * constructs a RDPScoring object with the same parameters as the given one
	 * @param arg the RDPScore which parameters shall be used
	 */
	public RDPScoring(RDPScoring arg) {
		this(
				arg.gamma, arg.delta, arg.epsilon, arg.zeta,
				arg.mutationMatrix,null
		);
	}
	
	
	public void setVorobin(String vorobin) {
		this.vorobin = vorobin;
	}

	/**
	 * calculates the score for a given OR node
	 * \phi (f, A, B) = \gamma * \phi^S(f,A,B) +	// mutation matrix (e.g. DAYHOFF)
	 * 					\delta * \phi^C(f,A,B) +	// contact capacity potential (see 123D)
	 * 					\epsilon * \phi^H(f,A,B) +	// hydrophobicity
	 * 					\zeta * \phi^P(f,A,B) -		// pair interaction
	 * 					GAP(f,A,B)					// insertions and deletions
	 * @param node the OR node that must be scored
	 * @return the score for the OR node (or rather the node's alignment)
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
		
		double result = 0.0;
		
		char[][] rows = f.getRows();
		
		for (int pos = 0; pos < rows[0].length; pos++) {
			// if positions are aligned
			if((rows[0][pos] != '-') && (rows[1][pos] != '-')) {
				// sum up score from mutation matrix
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

		double result = 0.0;
		
		char[][] rows = f.getRows();
		
		for (int pos = 0; pos < rows[0].length; pos++) {
			// TODO
			
		}
		
		return result;
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
		
		double result = 0.0;
		
		char[][] rows = f.getRows();
		
		int temppos = 0;	// position in template
		int targpos = 0;	// position in target
		
		for (int pos = 0; pos < rows[0].length; pos++) {
			
			if (rows[0][pos] == '-') {
				targpos++;
			} else if ((rows[1][pos] != '-')) {
				temppos++;
			} else {
				// sum up score
				result += hydrophobicity(a.getComp(targpos)) * 
						dob(b, temppos);
				temppos++;
				targpos++;
			}
		}
		
		return result;
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
		
		double result = 0.0;
		
		char[][] rows = f.getRows();
		
		int temppos = 0;	// position in template
		int targpos = 0;	// position in target
		
		for (int pos = 0; pos < rows[0].length; pos++) {
			
			if (rows[0][pos] == '-') {
				targpos++;
				if (gap.aligned) {
					// TODO find which gap is a real gap (i.e. don't calculate on not-yet
					// aligned parts)
					// insertion: target has aa, template not.
					// TODO result += phiP(f, a, b) at position pos
				}
			} else if ((rows[1][pos] != '-')) {
				temppos++;
				if (gap.aligned) {
					// TODO find which gap is a real gap (i.e. don't calculate on not-yet
					// aligned parts)
					// insertion: target has aa, template not.
					// TODO result += phiP(f, a, b) at position pos
				}
				// TODO result += phiC(f, a, b) at position pos
			} else {
				// result must not be changed here
				temppos++;
				targpos++;
			}
		}
		return result;
	}
	
	/**
	 * returns the hydrophobicity of an amino acid with given one letter code
	 * @param c
	 * @return the hydrophobicity of given amino acid
	 */
	private double hydrophobicity (char c) {
		// TODO
		
		return 0.0;
	}
	
	/**
	 * calculates the degree of burial (dob) for the given amino acid in the
	 * given structure
	 * @param structure the 3d structure of the template
	 * @param pos the position of the amino acid in the template
	 * @return the degree of burial [0..1]
	 */
	private double dob(PDBEntry structure, int pos) {
		return dob(structure, pos, 8.9, 1.0, 6.5);
	}
	
	/**
	 * calculates the degree of burial (dob) for the given amino acid in the
	 * given structure
	 * @param structure the 3d structure of the template
	 * @param pos the position of the amino acid in the template
	 * @param gridExtend value in Angstrom, additional space which will be filled by grid, CAVE: MUST be greater then gridClash!!
	 * @param gridDensisty value in Angstrom, denisty of solvent points with which grid will be filled
	 * @param gridClash value in Angstrom, distance every solvent must have to every peptide atom!
	 * @return the degree of burial [0..1]
	 */
	private double dob(PDBEntry structure, int pos, double gridExtend, double gridDensity, double gridClash) {
		VoroPPWrap voro = new VoroPPWrap(vorobin);
		VoronoiData data = new VoronoiData(structure.getID());
		data.reducePDB(VoroPrepType.CA, structure);
		data.fillGridWithoutClashes(8.9, 1.0, 6.5);
		voro.decomposite(data);
		HashMap<Integer,Double> faces = data.getFaces().get(pos);
		Set<Integer> solvents = data.getOuterGridIds();
		double outer = 0.0d;
		double inner = 0.0d;
		for(int neighbor : faces.keySet()){
			if(solvents.contains(neighbor)){
				outer += faces.get(neighbor);
			}else{
				inner += faces.get(neighbor);
			}
		}
		return outer/(outer+inner);
	}

}

/******************************************************************************
 * "A question that sometimes drives me hazy:                                 *
 *  Am I or are the others crazy?"                                            *
 *     - Albert Einstein (1979 - 1966)                                        *
 ******************************************************************************/
