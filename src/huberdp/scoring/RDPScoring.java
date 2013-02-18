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
import bioinfo.energy.potential.hydrophobicity.HydrophobicityMatrix;
import bioinfo.energy.potential.voronoi.VoroPPWrap;
import bioinfo.energy.potential.voronoi.VoroPrepType;
import bioinfo.energy.potential.voronoi.VoronoiData;
import bioinfo.proteins.PDBEntry;
import huberdp.RDPSolutionTreeOrNode;
import huberdp.RDPSolutionTreeAndNode;
import huberdp.Scoring;

/**
 * RDPScoring is an implementation of the scoring function given in the paper.
 * (Protein Threading by Recursive Dynamic Programming. JMB 290, 757-779)
 * 
 * @author huberste
 * @lastchange 2013-02-18
 */
public class RDPScoring implements Scoring {

	// TODO optimize scoring. Only one for() is needed.

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
	 * static reference to voro++ path
	 */
	private final static String VOROPATH = "/home/h/huberste/gobi/tina/tools/voro++_ubuntuquantal";

	/**
	 * empirically calibratet value for voro++
	 */
	private final static double GRID_EXTEND = 8.9;
	/**
	 * empirically calibratet value for voro++
	 */
	private final static double GRID_DENSITY = 1.0;
	/**
	 * empirically calibratet value for voro++
	 */
	private final static double GRID_CLASH = 6.5;
	/**
	 * empirically calibrate value for voro++
	 */
	private final static double MIN_CONTACT = 2.0;

	/**
	 * weight of the mutation matrix score
	 */
	private double gamma;
	/**
	 * weight of the contact capacity score
	 */
	private double delta;
	/**
	 * weight of the hydrophobicity score
	 */
	private double epsilon;
	/**
	 * weight of the pair interaction score
	 */
	private double zeta;

	/**
	 * mutation matrix for phiS
	 */
	private double[][] mutationMatrix;

	/**
	 * for the hydrophobicity scoring part phiH
	 */
	private HydrophobicityMatrix hydrophobicityMatrix;

	/**
	 * structure of the template
	 */
	private PDBEntry templateStructure;

	/**
	 * absolute path to location of voro++ binary (dont have one? look at
	 * ./tools/voro++)
	 * 
	 * @see http://math.lbl.gov/voro++/download/
	 */
	private String vorobin;

	/**
	 * Voronoi constants
	 */
	private double gridExtend;
	private double gridDensity;
	private double gridClash;
	private double minContact;

	/**
	 * Voro++ stuff
	 */
	private VoroPPWrap voro;
	private VoronoiData data;
	private Set<Integer> solvents;

	/**
	 * constructs a RDPScoring object with given paramters
	 * 
	 * @param gamma
	 *            weight of the mutation matrix score
	 * @param delta
	 *            weight of the contact capacity score
	 * @param epsilon
	 *            weight of the hydrophobicity score
	 * @param zeta
	 *            weight of the pair interaction score
	 * @param mutationMatrix
	 *            the mutation matrix that is to be used
	 * @param templatestructure
	 *            the template's structure
	 * @param vorobin
	 *            absolute path to location of voro++ binary
	 * @param gridExtend
	 *            variable vor voro++
	 * @param gridDensity
	 *            variable vor voro++
	 * @param gridClash
	 *            variable vor voro++
	 * @param minContact
	 *            variable vor voro++
	 */
	public RDPScoring(double gamma, double delta, double epsilon, double zeta,
			double[][] mutationMatrix,
			HydrophobicityMatrix hydrophobicityMatrix,
			PDBEntry templatestructure, String vorobin, double gridExtend,
			double gridDensity, double gridClash, double minContact) {
		this.gamma = gamma;
		this.delta = delta;
		this.epsilon = epsilon;
		this.zeta = zeta;
		this.mutationMatrix = mutationMatrix;
		this.hydrophobicityMatrix = hydrophobicityMatrix;
		this.vorobin = vorobin;
		setVoroVars(gridExtend, gridDensity, gridClash, minContact);
	}

	/**
	 * construcs a RDPScoring object with standard parameters
	 */
	public RDPScoring() {
		this(GAMMA, DELTA, EPSILON, ZETA,
				bioinfo.alignment.matrices.QuasarMatrix.DAYHOFF_MATRIX,
				new HydrophobicityMatrix(), null, VOROPATH, GRID_EXTEND,
				GRID_DENSITY, GRID_CLASH, MIN_CONTACT);
	}

	/**
	 * constructs a RDPScoring object with the same parameters as the given one
	 * 
	 * @param arg
	 *            the RDPScore which parameters shall be used
	 */
	public RDPScoring(RDPScoring arg) {
		this(arg.gamma, arg.delta, arg.epsilon, arg.zeta, arg.mutationMatrix,
				arg.hydrophobicityMatrix, arg.templateStructure, arg.vorobin,
				arg.gridExtend, arg.gridDensity, arg.gridClash, arg.minContact);
	}

	/**
	 * initializes Voro++ stuff
	 * 
	 * @param gridExtend
	 *            value in Angstrom, additional space which will be filled by
	 *            grid, CAVE: MUST be greater then gridClash!!
	 * @param gridDensisty
	 *            value in Angstrom, denisty of solvent points with which grid
	 *            will be filled
	 * @param gridClash
	 *            value in Angstrom, distance every solvent must have to every
	 *            peptide atom!
	 */
	public void initVoro() {
		if (vorobin != null && templateStructure != null) {
			voro = new VoroPPWrap(vorobin);
			data = new VoronoiData(templateStructure.getID());
			data.reducePDB(VoroPrepType.CA, templateStructure);
			data.fillGridWithoutClashes(gridExtend, gridDensity, gridClash);
			voro.decomposite(data);
			data.detectOuterGrid(minContact);
			solvents = data.getOuterGridIds();
		}
	}

	/**
	 * sets binary path of voro++
	 * 
	 * @param vorobin
	 * @author seitza
	 */
	public void setVorobin(String vorobin) {
		this.vorobin = vorobin;
	}

	/**
	 * 
	 * @param gridExtend
	 * @param gridDensity
	 * @param gridClash
	 * @param minContact
	 */
	public void setVoroVars(double gridExtend, double gridDensity,
			double gridClash, double minContact) {
		this.gridExtend = gridExtend;
		this.gridDensity = gridDensity;
		this.gridClash = gridClash;
		this.minContact = minContact;
	}

	/**
	 * calculates the score for a given OR node \phi (f, A, B) = \gamma *
	 * \phi^S(f,A,B) + // mutation matrix (e.g. DAYHOFF) \delta * \phi^C(f,A,B)
	 * + // contact capacity potential (see 123D) \epsilon * \phi^H(f,A,B) + //
	 * hydrophobicity \zeta * \phi^P(f,A,B) - // pair interaction GAP(f,A,B) //
	 * insertions and deletions
	 * 
	 * @param node
	 *            the OR node that must be scored
	 * @return the score for the OR node (or rather the node's alignment)
	 */
	@Override
	public double score(RDPSolutionTreeOrNode node) {

		double result = 0.0;

		// check if correct structure is set
		if ((this.templateStructure == null)
				|| (templateStructure != node.getProblem().templateStructure)) {
			templateStructure = node.getProblem().templateStructure;
			initVoro();
		}

		// check if voronoi composition is set
		/*
		 * // normally this should never be the case. if (voro == null) {
		 * initVoro(); }
		 */

		if (node.getParent() != null) { // node is not root
			// add parent's alignment's score to parent's parent's score
			result += ((RDPSolutionTreeOrNode) node.getParent().getParent())
					.getScore();

			SequenceAlignment f = ((RDPSolutionTreeAndNode) node.getParent())
					.getPA().alignment;
			Sequence a = node.getProblem().targetSequence;
			PDBEntry b = node.getProblem().templateStructure;

			result = gamma * phiS(f, a, b) + delta * phiC(f, a, b) + epsilon
					* phiH(f, a, b) + zeta * phiP(f, a, b) - gap(f, a, b);
		}

		return result;
	}

	/**
	 * "phiS scores the alignment f with respect to well-known sequence based
	 * mutation matrices" (From: Protein Threading by Recursive Dynamic
	 * Programming. JMB 290, 757-779)
	 * 
	 * @param f
	 *            the alignment (so far)
	 * @param a
	 *            the target sequence
	 * @param b
	 *            the template structure
	 * @return the calculated sequence-based score
	 */
	private double phiS(SequenceAlignment f, Sequence a, PDBEntry b) {

		double result = 0.0;

		char[][] rows = f.getRows();

		for (int pos = 0; pos < rows[0].length; pos++) {
			// if positions are aligned
			if ((rows[0][pos] != '-') && (rows[1][pos] != '-')) {
				// sum up score from mutation matrix
				result += mutationMatrix[rows[0][pos] - 65][rows[1][pos] - 65];
			}
		}

		return result;
	}

	/**
	 * "(...) contact-capacity-potential phiC (...)" (From: Protein Threading by
	 * Recursive Dynamic Programming. JMB 290, 757-779)
	 * 
	 * @param f
	 *            the alignment (so far)
	 * @param a
	 *            the target sequence
	 * @param b
	 *            the template structure
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
	 * "phiH [scores] (...) the hydrophobicity (...)" (From: Protein Threading
	 * by Recursive Dynamic Programming. JMB 290, 757-779)
	 * 
	 * @param f
	 *            the alignment (so far)
	 * @param a
	 *            the target sequence
	 * @param b
	 *            the template structure
	 * @return the calculated hydrophobicity based score
	 */
	private double phiH(SequenceAlignment f, Sequence a, PDBEntry b) {

		double result = 0.0;

		char[][] rows = f.getRows();

		int buckets = hydrophobicityMatrix.getBuckets();

		int temppos = 0; // position in template
		int targpos = 0; // position in target

		for (int pos = 0; pos < rows[0].length; pos++) {

			if (rows[0][pos] == '-') { // if insertion
				targpos++;
			} else if ((rows[1][pos] == '-')) { // if deletion
				temppos++;
			} else { // if match
				// sum up score
				int astype = (a.getComp(targpos)) - 65;
				double dob = dob(b, temppos);
				for (int bucket = 0; bucket < buckets; bucket++) {
					if (dob <= ((double) bucket + 1.0)
							* (1.0 / (double) buckets)) {
						result += hydrophobicityMatrix.getValue(astype, bucket);
						break;
					}
				}
				temppos++;
				targpos++;
			}
		}

		return result;
	}

	/**
	 * "phiP denotes the pair interaction term (...)" (From: Protein Threading
	 * by Recursive Dynamic Programming. JMB 290, 757-779)
	 * 
	 * @param f
	 *            the alignment (so far)
	 * @param a
	 *            the target sequence
	 * @param b
	 *            the template structure
	 * @return the calculated pair interaction based score
	 */
	private double phiP(SequenceAlignment f, Sequence a, PDBEntry b) {
		// TODO
		return 0.0;
	}

	/**
	 * "GAP penalizes insertions and deletions." (From: Protein Threading by
	 * Recursive Dynamic Programming. JMB 290, 757-779)
	 * 
	 * @param f
	 *            the alignment (so far)
	 * @param a
	 *            the target sequence
	 * @param b
	 *            the template structure
	 * @return the calculated pair interaction based score
	 */
	private double gap(SequenceAlignment f, Sequence a, PDBEntry b) {

		double result = 0.0;

		char[][] rows = f.getRows();

		int temppos = 0; // position in template
		int targpos = 0; // position in target

		for (int pos = 0; pos < rows[0].length; pos++) {

			if (rows[0][pos] == '-') {
				targpos++;
				if (true) {
					// TODO find if gap is a real gap (i.e. don't calculate on
					// not-yet
					// aligned parts)
					// insertion: target has aa, template not.
					// TODO result += phiP(f, a, b) at position pos
				}
			} else if ((rows[1][pos] != '-')) {
				temppos++;
				if (true) {
					// TODO find if gap is a real gap (i.e. don't calculate on
					// not-yet
					// aligned parts)
					// deletion: template has aa, target not.
					// TODO result += phiC(f, a, b) at position pos
				}
				//
			} else {
				// result need not to be changed here
				temppos++;
				targpos++;
			}
		}
		return result;
	}

	/**
	 * calculates the degree of burial (dob) for the given amino acid in the
	 * given structure
	 * 
	 * @author seitza
	 * @param structure
	 *            the 3d structure of the template
	 * @param pos
	 *            the position of the amino acid in the template
	 * @return the degree of burial [0..1]
	 */
	private double dob(PDBEntry structure, int pos) {

		double outer = 0.0;
		double inner = 0.0;

		HashMap<Integer, Double> faces = data.getFaces().get(pos);

		for (int neighbor : faces.keySet()) {
			if (solvents.contains(neighbor)) {
				outer += faces.get(neighbor);
			} else {
				inner += faces.get(neighbor);
			}
		}
		return outer / (outer + inner);
	}

}

/******************************************************************************
 * "A question that sometimes drives me hazy: * Am I or are the others crazy?" *
 * - Albert Einstein (1879 - 1955) *
 ******************************************************************************/
