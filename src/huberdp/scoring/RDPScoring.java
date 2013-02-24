/******************************************************************************
 * huberdp.scoring.RDPScoring.java                                            *
 *                                                                            *
 * This file contains the class RDPScoring which is RDP's scoring function.   *
 *                                                                            *
 * This file is best read at line width 80 and tab width 4.                   *
 *                                                                   huberste *
 ******************************************************************************/
package huberdp.scoring;

import static util.Util.*;

import java.util.HashMap;
import java.util.Set;

import bioinfo.Sequence;
import bioinfo.alignment.Threading;
import bioinfo.energy.potential.PotentialDimension;
import bioinfo.energy.potential.SipplContactPotential;
import bioinfo.energy.potential.hydrophobicity.HydrophobicityMatrix;
import bioinfo.energy.potential.voronoi.VoroPPWrap;
import bioinfo.energy.potential.voronoi.VoroPrepType;
import bioinfo.energy.potential.voronoi.VoronoiData;
import bioinfo.proteins.AminoAcid;
import bioinfo.proteins.CCPMatrix;
import bioinfo.proteins.DSSPEntry;
import bioinfo.proteins.DSSPFileReader;
import bioinfo.proteins.PDBEntry;
import bioinfo.proteins.SecStructEight;
import bioinfo.proteins.SecStructThree;
import huberdp.Scoring;

/**
 * RDPScoring is an implementation of the scoring function given in the paper.
 * (Protein Threading by Recursive Dynamic Programming. JMB 290, 757-779)
 * 
 * @author huberste
 * @lastchange 2013-02-22
 */
public class RDPScoring implements Scoring {

	/**
	 * empirically calibrated weights for the scoring function
	 */
	public final static double GAMMA = 1.0, DELTA = 0.1, EPSILON = 2.0,
			ZETA = 4.0, GAP = 14.0;

	/**
	 * empirically calibrated values for voro++
	 */
	public final static double GRID_EXTEND = 8.9, GRID_DENSITY = 1.0,
			GRID_CLASH = 6.5, MIN_CONTACT = 2.0;

	/**
	 * weights for the scoring function
	 */
	private double gamma, delta, epsilon, zeta, gap_weight;

	/**
	 * mutation matrix for phiS
	 */
	private double[][] mutationMatrix;

	/**
	 * for the hydrophobicity scoring part phiH
	 */
	private HydrophobicityMatrix hydrophobicityMatrix;

	/**
	 * ContactCapacityMatrix
	 */
	private CCPMatrix ccpMatrix;

	/**
	 * Folder- and filenames
	 */
	private String dsspFolder, vorobin, tempdir;

	/**
	 * path to the vpot file (PairPotential File)
	 */
	private SipplContactPotential pcp;

	/**
	 * structure of the template
	 */
	private PDBEntry templateStructure;

	/**
	 * Voronoi constants
	 */
	private double gridExtend, gridDensity, gridClash, minContact;

	/**
	 * Voro++ stuff
	 */
	private VoroPPWrap voro;
	/**
	 * Voro++ stuff
	 */
	private VoronoiData data;
	/**
	 * Voro++ stuff
	 */
	private Set<Integer> solvents;

	private int[][] contacts;

	private SecStructThree[] ss;

	private double[] avgPotential;

	/**
	 * for HydrobhobicityMatrix
	 */
	private int hydroBuckets;

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
	 * @param gap
	 *            weight of the gap score
	 * @param mutationMatrix
	 *            the mutation matrix that is to be used
	 * @param hydrophobicityMatrix
	 *            the hydrophobicityMatrix to be used
	 * @param ccpMatrix
	 *            the CCPMatrix to be used
	 * @param dsspFolder
	 *            folder to the dssp files
	 * @param sippl
	 *            sippl Pair Potential
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
	 * @param tempdir
	 *            temporary directory for vorolign
	 */
	public RDPScoring(double gamma, double delta, double epsilon, double zeta,
			double gap, double[][] mutationMatrix,
			HydrophobicityMatrix hydrophobicityMatrix, CCPMatrix ccpMatrix,
			String dsspFolder, SipplContactPotential sippl,
			PDBEntry templatestructure, String vorobin, double gridExtend,
			double gridDensity, double gridClash, double minContact,
			String tempdir) {
		this.gamma = gamma;
		this.delta = delta;
		this.epsilon = epsilon;
		this.zeta = zeta;
		this.gap_weight = -gap;
		this.mutationMatrix = mutationMatrix;
		this.hydrophobicityMatrix = hydrophobicityMatrix;
		this.pcp = sippl;
		this.ccpMatrix = ccpMatrix;
		this.dsspFolder = dsspFolder;
		this.vorobin = vorobin;
		this.gridExtend = gridExtend;
		this.gridDensity = gridDensity;
		this.gridClash = gridClash;
		this.minContact = minContact;
		this.tempdir = tempdir;
		calculateAveragePPs();
	}

	/**
	 * constructs a RDPScoring object with the same parameters as the given one
	 * 
	 * @param arg
	 *            the RDPScore which parameters shall be used
	 */
	public RDPScoring(RDPScoring arg) {
		this(arg.gamma, arg.delta, arg.epsilon, arg.zeta, arg.gap_weight,
				arg.mutationMatrix, arg.hydrophobicityMatrix, arg.ccpMatrix,
				arg.dsspFolder, arg.pcp, arg.templateStructure, arg.vorobin,
				arg.gridExtend, arg.gridDensity, arg.gridClash, arg.minContact,
				arg.tempdir);
	}

	/**
	 * calculates AveragePPs for insertions
	 */
	private void calculateAveragePPs() {
		avgPotential = new double[26];
		PotentialDimension<Double> potential = pcp.getPotential();
		for (int aa1 = 0; aa1 < 26; aa1++) {
			for (int k = 0; k < 6; k++) {
				for (int d = 0; d < 6; d++) {
					for (int aa2 = 0; aa2 < 26; aa2++) {
						int[] path1 = { k, d, aa1, aa2 };
						avgPotential[aa1] += potential.getByAddress(path1)
								.getValue();
						int[] path2 = { k, d, aa2, aa1 };
						avgPotential[aa1] += potential.getByAddress(path2)
								.getValue();
					}
				}
			}
			avgPotential[aa1] = avgPotential[aa1] / 936.0; // 936 = 26*6*6
		}
	}

	/**
	 * initializes stuff
	 */
	public void init() {
		// set & calculate Voro++ stuff
		if (vorobin != null && templateStructure != null) {
			voro = new VoroPPWrap(tempdir, vorobin);
			data = new VoronoiData(templateStructure.getID());
			data.reducePDB(VoroPrepType.CA, templateStructure);
			data.fillGridWithoutClashes(gridExtend, gridDensity, gridClash);
			voro.decomposite(data);
			data.detectOuterGrid(minContact);
			solvents = data.getOuterGridIds();
		}

		// calculate contact matrix
		if (templateStructure != null) {
			// Count contacts
			// first dimension: local / long range
			// second dimension: position in structure
			contacts = new int[2][templateStructure.length()];
			// calculate contacts fo every amino acid
			for (int partnera = 0; partnera < templateStructure.length(); partnera++) {
				for (int partnerb = partnera + 1; partnerb < templateStructure
						.length(); partnerb++) {
					if (calcDistance(templateStructure.getAminoAcid(partnera),
							templateStructure.getAminoAcid(partnerb)) < 7.0) {
						if (Math.abs(partnera - partnerb) < 5) { // local
							contacts[0][partnera]++;
							contacts[0][partnerb]++;
						} else { // longRange
							contacts[1][partnera]++;
							contacts[1][partnerb]++;
						}
					}
				}
			}

			// read SecStruct from DSSP File
			String dsspFileName = dsspFolder
					+ templateStructure.getId().toLowerCase()
					+ templateStructure.getID().substring(4, 7) + ".dssp";
			DSSPEntry dssp = DSSPFileReader.readDSSPFile(dsspFileName);
			SecStructEight[] temp = dssp.getSecondaryStructure();
			ss = new SecStructThree[temp.length];
			for (int i = 0; i < temp.length; i++) {
				ss[i] = temp[i].getThreeClassAnalogon();
			}
		}

		// initialize other variables
		hydroBuckets = hydrophobicityMatrix.getBuckets();

	}

	/**
	 * calculates the score for a given Threading \phi (f, A, B) = \gamma *
	 * \phi^S(f,A,B) + // mutation matrix (e.g. DAYHOFF) \delta * \phi^C(f,A,B)
	 * + // contact capacity potential (see 123D) \epsilon * \phi^H(f,A,B) + //
	 * hydrophobicity \zeta * \phi^P(f,A,B) - // pair interaction GAP(f,A,B) //
	 * insertions and deletions
	 * 
	 * @param threading
	 *            the Threading to be scored
	 * @return the score for the Threading
	 */
	@Override
	public double score(Threading threading) {

		// check if correct structure is set
		if ((this.templateStructure == null)
				|| (templateStructure != threading.getStructure())) {
			templateStructure = threading.getStructure();
			init();
		}

		double sscore = 0.0;
		double cscore = 0.0;
		double hscore = 0.0;
		double pscore = 0.0;
		double gapscore = 0.0;

		char[][] rows = threading.getRowsAsCharArray();
		int temppos = 0; // position in template
		int targpos = 0; // position in target

		PDBEntry b = threading.getStructure();

		for (int pos = 0; pos < rows[0].length; pos++) {

			if (rows[0][pos] == '-') { // if insertion
				// TODO find if gap is a real gap (i.e. don't calculate on
				// not-yet aligned parts)
				gapscore += getInsertionScore(threading, targpos);
				// count up
				targpos++;
			} else if ((rows[1][pos] == '-')) { // if deletion
				// TODO find if gap is a real gap (i.e. don't calculate on
				// not-yet aligned parts)
				gapscore += getDeletionScore(threading, temppos);
				// count up
				temppos++;
			} else { // if match
				int astemp = rows[0][pos] - 65;
				int astarg = rows[1][pos] - 65;
				// phiS
				sscore += mutationMatrix[astemp][astarg];
				// phiC
				cscore += ccpMatrix.getValue(astarg, ss[temppos], 0,
						contacts[0][temppos]);
				cscore += ccpMatrix.getValue(astarg, ss[temppos], 1,
						contacts[1][temppos]);
				// phiH
				double dob = dob(b, temppos);
				for (int bucket = 0; bucket < hydroBuckets; bucket++) {
					if (dob <= ((double) bucket + 1.0) / (double) hydroBuckets) {
						hscore += hydrophobicityMatrix.getValue(astarg, bucket);
						break;
					}
				}
				// count up
				temppos++;
				targpos++;
			}
		}

		// phiP
		PDBEntry model = modifyModel(threading);
		pscore = pcp.scoreModel(model);

		// hint: gap_weight is applicated in getInsertionScore /
		// getDeletionScore
		return gamma * sscore + delta * cscore + epsilon * hscore + zeta
				* pscore + gapscore;
	}

	/**
	 * 
	 * @param threading
	 *            the treading to be scored
	 * @param m
	 *            the position in the template
	 * @param n
	 *            the position in the target
	 * @return the score for the aligned positions n (target) and m (template)
	 */
	public double getScore(Threading threading, int m, int n) {

		// check if correct structure is set
		if ((this.templateStructure == null)
				|| (templateStructure != threading.getStructure())) {
			templateStructure = threading.getStructure();
			init();
		}

		return gamma * phiS(threading, m, n) + delta * phiC(threading, m, n)
				+ epsilon * phiH(threading, m, n) + zeta
				* phiP(threading, m, n);

	}

	/**
	 * "phiS scores the alignment f with respect to well-known sequence based
	 * mutation matrices" (From: Protein Threading by Recursive Dynamic
	 * Programming. JMB 290, 757-779)
	 * 
	 * @param f
	 *            the alignment (so far)
	 * @param m
	 *            position in the template
	 * @param n
	 *            position in th target
	 * @return the calculated sequence-based score
	 */
	private double phiS(Threading f, int m, int n) {

		int a = f.getStructure().getAminoAcid(m).getName().getNumber();
		int b = f.getSequence().getComp(n) - 65;

		return mutationMatrix[a][b];
	}

	/**
	 * "(...) contact-capacity-potential phiC (...)" (From: Protein Threading by
	 * Recursive Dynamic Programming. JMB 290, 757-779) <br />
	 * Two AminoAcids are in contact if their C alpha atoms are less than 7 Å
	 * distant
	 * 
	 * @param f
	 *            the alignment (so far)
	 * @param m
	 *            position in the template
	 * @param n
	 *            position in th target
	 * @return the calculated contact-capacity based score
	 */
	private double phiC(Threading f, int m, int n) {

		int aa = f.getSequence().getComp(n) - 65;

		double result = ccpMatrix.getValue(aa, ss[m], 0, contacts[0][m])
				+ ccpMatrix.getValue(aa, ss[m], 1, contacts[1][m]);

		return result;
	}

	/**
	 * "phiH [scores] (...) the hydrophobicity (...)" (From: Protein Threading
	 * by Recursive Dynamic Programming. JMB 290, 757-779)
	 * 
	 * @param f
	 *            the threading
	 * @param m
	 *            the position in the template
	 * @param n
	 *            the position in the target
	 * @return the calculated hydrophobicity based score
	 */
	private double phiH(Threading f, int m, int n) {

		int astype = (f.getSequence().getComp(n)) - 65;
		double dob = dob(f.getStructure(), m);
		for (int bucket = 0; bucket < hydroBuckets; bucket++) {
			if (dob <= ((double) bucket + 1.0) / (double) hydroBuckets) {
				return hydrophobicityMatrix.getValue(astype, bucket);
			}
		}

		return 0.0;
	}

	/**
	 * "phiP denotes the pair interaction term (...)" (From: Protein Threading
	 * by Recursive Dynamic Programming. JMB 290, 757-779)
	 * 
	 * @param f
	 *            the alignment (so far)
	 * @param m
	 *            position in the template
	 * @param n
	 *            position in the target
	 * @return the calculated pair interaction based score
	 */
	private double phiP(Threading f, int m, int n) {
		// use SipplContactPotential from bioinfo.energy.potential
		PDBEntry model = modifyModel(f);
		return pcp.getAminoScores(model)[m];
	}

	/**
	 * calculates the degree of burial (dob) for the given amino acid in the
	 * given structure
	 * 
	 * @author seitza
	 * @author huberste
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

	/**
	 * modifies the model so it can be scored by phiP. It replaces all aligned
	 * AminoAcids.
	 * 
	 * @param f
	 *            SequenceAlignment (template, target)
	 * @return a new, modified model
	 */
	private static PDBEntry modifyModel(Threading f) {

		char[][] rows = f.getRowsAsCharArray();
		PDBEntry b = f.getStructure();
		Sequence a = f.getSequence();

		AminoAcid[] aminoAcids = new AminoAcid[b.length()];

		int temppos = 0;
		int targpos = 0;
		for (int pos = 0; pos < rows[0].length; pos++) {

			if (rows[0][pos] == '-') { // if insertion
				targpos++;
			} else if ((rows[1][pos] == '-')) { // if deletion
				aminoAcids[temppos] = b.getAminoAcid(temppos);
				temppos++;
			} else { // if match
				aminoAcids[temppos] = new AminoAcid(a.getComp(targpos), b
						.getAminoAcid(temppos).getResIndex(), b.getAminoAcid(
						temppos).getAtoms());
				temppos++;
				targpos++;
			}
		}

		PDBEntry result = new PDBEntry(b.getId(), b.getChainID(),
				b.getChainIDNum(), aminoAcids);

		return result;
	}

	/**
	 * Calculates the score for an insertion of target's AminoAcid n.
	 * 
	 * @param threading
	 *            the threading to be scored
	 * @param n
	 *            the number of the AminoAcid in the target that is inserted
	 * @return the average PairPotential of the inserted AA
	 */
	public double getInsertionScore(Threading threading, int n) {
		// hint: gap_weight gets negative in constructor!
		return gap_weight
				* (avgPotential[threading.getSequence().getComp(n) - 65]);
	}

	/**
	 * Calculates the score for an deletion of template's AminoAcid m.
	 * 
	 * @param threading
	 *            the threading to be scored
	 * @param m
	 *            the number of the AminoAcid in the template that is deleted
	 * @return the average PairPotential of the deleted AA
	 */
	public double getDeletionScore(Threading threading, int m) {
		// hint: gap_weight gets negavite in constructor!
		return gap_weight
				* (avgPotential[threading.getStructure().getAminoAcid(m)
						.getName().getNumber()]);
	}

}

/******************************************************************************
 * "A question that sometimes drives me hazy:                                 *
 *  Am I or are the others crazy?"                                            *
 * - Albert Einstein (1879 - 1955)                                            *
 ******************************************************************************/
