/******************************************************************************
 * huberdp.Scoring.RDPScoring.java                                            *
 * This file contains the class RDPScoring which is RDP's scoring function.   *
 *                                                                            *
 * This file is best read at line width 80 and tab width 4.                   *
 *                                                                   huberste *
 ******************************************************************************/
package huberdp.scoring;

import huberdp.Scoring;

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
import bioinfo.proteins.Atom;
import bioinfo.proteins.AtomType;
import bioinfo.proteins.CCPMatrix;
import bioinfo.proteins.DSSPEntry;
import bioinfo.proteins.DSSPFileReader;
import bioinfo.proteins.PDBEntry;
import bioinfo.proteins.SecStructEight;
import bioinfo.proteins.SecStructThree;

/**
 * RDPScoring is an implementation of the scoring function given in the paper.
 * (Protein Threading by Recursive Dynamic Programming. JMB 290, 757-779)
 * 
 * @author huberste
 * @lastchange 2013-02-18
 */
public class RDPScoring implements Scoring {

	public double smax, smin, cmax, cmin, hmax, hmin, pmax, pmin, scoremin,
			scoremax, deletemin, deletemax, insertmin, insertmax;
	public double ssum, csum, hsum, psum, scoresum, deletesum, insertsum;
	public int num, deletenum, insertnum;

	// TODO optimize scoring. Only one for() is needed.

	// TODO calibrate Scoring weights

	/**
	 * empirically calibrated weight of the mutation matrix score
	 */
	public final static double GAMMA = 1.0;
	/**
	 * empirically calibrated weight of the contact capacity score
	 */
	public final static double DELTA = 0.1;
	/**
	 * empirically calibrated weight of the hydrophobicity score
	 */
	public final static double EPSILON = 2.0;
	/**
	 * empirically calibrated weight of the pair interaction score
	 */
	public final static double ZETA = 4.0;
	/**
	 * empirically calibrated weight of the gaps
	 */
	public final static double GAP = 3.0;

	/**
	 * static reference to voro++ path
	 */
	public final static String VOROPATH = "/home/h/huberste/gobi/tina/tools/voro++_ubuntuquantal";

	/**
	 * empirically calibratet value for voro++
	 */
	public final static double GRID_EXTEND = 8.9;
	/**
	 * empirically calibratet value for voro++
	 */
	public final static double GRID_DENSITY = 1.0;
	/**
	 * empirically calibratet value for voro++
	 */
	public final static double GRID_CLASH = 6.5;
	/**
	 * empirically calibrate value for voro++
	 */
	public final static double MIN_CONTACT = 2.0;

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
	 * additionally weight parameter for GAPs
	 */
	private double gap_weight;

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

	private String dsspFolder;

	/**
	 * path to the vpot file (PairPotential File)
	 */
	private SipplContactPotential pcp;

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
	 * Voronoi constant
	 */
	private double gridExtend;
	/**
	 * Voronoi constant
	 */
	private double gridDensity;
	/**
	 * Voronoi constant
	 */
	private double gridClash;
	/**
	 * Voronoi constant
	 */
	private double minContact;

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
	 * @param hydrophobicityMatrix
	 *            the hydrophobicityMatrix to be used
	 * @param ccpMatrix
	 *            the CCPMatrix to be used
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
			double gap, double[][] mutationMatrix,
			HydrophobicityMatrix hydrophobicityMatrix, CCPMatrix ccpMatrix,
			String dsspFolder, SipplContactPotential sippl,
			PDBEntry templatestructure, String vorobin, double gridExtend,
			double gridDensity, double gridClash, double minContact) {
		this.gamma = gamma;
		this.delta = delta;
		this.epsilon = epsilon;
		this.zeta = zeta;
		this.gap_weight = gap;
		this.mutationMatrix = mutationMatrix;
		this.hydrophobicityMatrix = hydrophobicityMatrix;
		this.pcp = sippl;
		this.ccpMatrix = ccpMatrix;
		this.dsspFolder = dsspFolder;
		this.vorobin = vorobin;
		calculateAveragePPs();
		setVoroVars(gridExtend, gridDensity, gridClash, minContact);
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
	public void init() {
		if (vorobin != null && templateStructure != null) {
			voro = new VoroPPWrap(vorobin);
			data = new VoronoiData(templateStructure.getLongID());
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
		}

		// read SecStruct from DSSP File
		String dsspFileName = dsspFolder
				+ templateStructure.getID().toLowerCase()
				+ templateStructure.getLongID().substring(4, 7) + ".dssp";
		DSSPEntry dssp = DSSPFileReader.readDSSPFile(dsspFileName);
		SecStructEight[] temp = dssp.getSecondaryStructure();
		ss = new SecStructThree[temp.length];
		for (int i = 0; i < temp.length; i++) {
			ss[i] = temp[i].getThreeClassAnalogon();
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
	 * sets voronoi variables
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

		double result = 0.0;

		// check if correct structure is set
		if ((this.templateStructure == null)
				|| (templateStructure != threading.getStructure())) {
			templateStructure = threading.getStructure();
			init();
		}

		result = gamma * phiS(threading) + delta * phiC(threading) + epsilon
				* phiH(threading) + zeta * phiP(threading) - gap(threading);

		return result;
	}

	/**
	 * "phiS scores the alignment f with respect to well-known sequence based
	 * mutation matrices" (From: Protein Threading by Recursive Dynamic
	 * Programming. JMB 290, 757-779)
	 * 
	 * @param f
	 *            the threading
	 * @return the calculated sequence-based score
	 */
	private double phiS(Threading f) {

		double result = 0.0;

		char[][] rows = f.getRowsAsCharArray();

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

		char a = f.getStructure().getAminoAcid(m).getName()
				.getThreeLetterCode().charAt(0);
		char b = f.getSequence().getComp(n);

		return mutationMatrix[a - 65][b - 65];
	}

	/**
	 * "(...) contact-capacity-potential phiC (...)" (From: Protein Threading by
	 * Recursive Dynamic Programming. JMB 290, 757-779) <br />
	 * Two AminoAcids are in contact if their C alpha atoms are less than 7 Å
	 * distant
	 * 
	 * @param f
	 *            the threading
	 * @return the calculated contact-capacity based score
	 */
	private double phiC(Threading f) {

		double result = 0.0;

		char[][] rows = f.getRowsAsCharArray();
		Sequence a = f.getSequence();

		// read SecStruct from DSSP File
		String dsspFileName = DSSPFileReader.DSSP_FOLDER
				+ f.getStructure().getID().toLowerCase()
				+ f.getStructure().getLongID().substring(4, 7) + ".dssp";
		DSSPEntry dssp = DSSPFileReader.readDSSPFile(dsspFileName);
		SecStructEight[] ss = dssp.getSecondaryStructure();

		int temppos = 0; // position in template
		int targpos = 0; // position in target

		for (int pos = 0; pos < rows[0].length; pos++) {

			if (rows[0][pos] == '-') { // if insertion
				targpos++;
			} else if ((rows[1][pos] == '-')) { // if deletion
				temppos++;
			} else { // if match
				// sum up score
				int tmp = contacts[0][temppos];
				result += ccpMatrix.getValue(a.getComp(targpos),
						ss[temppos].getThreeClassAnalogon(), 0, tmp);
				result += ccpMatrix.getValue(a.getComp(targpos),
						ss[temppos].getThreeClassAnalogon(), 1,
						contacts[1][temppos]);
				temppos++;
				targpos++;
			}
		}

		return result;
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

		char aa = f.getSequence().getComp(n);

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
	 * @return the calculated hydrophobicity based score
	 */
	private double phiH(Threading f) {

		double result = 0.0;

		char[][] rows = f.getRowsAsCharArray();
		PDBEntry b = f.getStructure();
		Sequence a = f.getSequence();

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

		double result = 0.0;

		PDBEntry b = f.getStructure();
		Sequence a = f.getSequence();

		int buckets = hydrophobicityMatrix.getBuckets();

		int astype = (a.getComp(n)) - 65;
		double dob = dob(b, m);
		for (int bucket = 0; bucket < buckets; bucket++) {
			if (dob <= ((double) bucket + 1.0) * (1.0 / (double) buckets)) {
				result = hydrophobicityMatrix.getValue(astype, bucket);
				break;
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
	private double phiP(Threading f) {
		// use SipplContactPotential from bioinfo.energy.potential
		PDBEntry model = modifyModel(f);
		return pcp.scoreModel(model);
		// return 0.0;
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
	private double gap(Threading f) {

		double result = 0.0;

		char[][] rows = f.getRowsAsCharArray();

		int temppos = 0; // position in template
		int targpos = 0; // position in target

		for (int pos = 0; pos < rows[0].length; pos++) {

			if (rows[0][pos] == '-') { // insertion: target has aa, template
										// not.
				// TODO find if gap is a real gap (i.e. don't calculate on
				// not-yet aligned parts)
				result += getInsertionScore(f, targpos);
				targpos++;
			} else if ((rows[1][pos] != '-')) { // deletion: template has aa,
												// target not.
				// TODO find if gap is a real gap (i.e. don't calculate on
				// not-yet aligned parts)
				result += getDeletionScore(f, temppos);
				temppos++;
			} else {
				// result need not to be changed here
				temppos++;
				targpos++;
			}
		}
		return result;
	}

	/**
	 * calculates the euklidian distance between two AminoAcids
	 * 
	 * @param a
	 *            an AmoniAcid
	 * @param b
	 *            another AminoAcid
	 * @return the euklidian distance between the two AminoAcid's C alpha atoms
	 */
	private double calcDistance(AminoAcid a, AminoAcid b) {
		Atom caa = a.getAtomByType(AtomType.CA);
		Atom cab = b.getAtomByType(AtomType.CA);
		if (caa != null && cab != null) {
			return calcDistance(caa, cab);
		}

		return 0.0;
	}

	/**
	 * calculates the distance between two atoms
	 * 
	 * @param a
	 *            an Atom
	 * @param b
	 *            another Atom
	 * @return the euklidian distance between two Atoms
	 */
	private double calcDistance(Atom a, Atom b) {
		double[] apos = a.getPosition();
		double[] bpos = b.getPosition();
		double[] dis = { apos[0] - bpos[0], apos[1] - bpos[1],
				apos[2] - bpos[2] };
		return Math.sqrt(Math.pow(dis[0], 2) + Math.pow(dis[1], 2)
				+ Math.pow(dis[2], 2));
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
	 * modifies the model so it can be scored by phiP
	 * 
	 * @param f
	 *            SequenceAlignment (template, target)
	 * @param a
	 *            target Sequence
	 * @param b
	 *            template structure
	 * @return
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

		PDBEntry result = new PDBEntry(b.getID(), b.getChainID(),
				b.getChainIDNum(), aminoAcids);

		return result;
	}

	/**
	 * 
	 * @param threading
	 * @param m
	 *            the position in the template
	 * @param n
	 *            the position in the target
	 * @return
	 */
	public double getScore(Threading threading, int m, int n) {

		// check if correct structure is set
		if ((this.templateStructure == null)
				|| (templateStructure != threading.getStructure())) {
			templateStructure = threading.getStructure();
			init();
		}

		double sscore = gamma * phiS(threading, m, n);
		smax = smax < sscore ? sscore : smax;
		smin = smin > sscore ? sscore : smin;
		ssum += sscore;
		double cscore = delta * phiC(threading, m, n);
		cmax = cmax < cscore ? cscore : cmax;
		cmin = cmin > cscore ? cscore : cmin;
		csum += cscore;
		double hscore = epsilon * phiH(threading, m, n);
		hmax = hmax < hscore ? hscore : hmax;
		hmin = hmin > hscore ? hscore : hmin;
		hsum += hscore;
		double pscore = zeta * phiP(threading, m, n);
		pmax = pmax < pscore ? pscore : pmax;
		pmin = pmin > pscore ? pscore : pmin;
		psum += pscore;

		double result = sscore + cscore + hscore + pscore;
		scoremax = scoremax < result ? result : scoremax;
		scoremin = scoremin > result ? result : scoremin;
		scoresum += result;

		num++;

		return result;

	}

	/**
	 * Calculates the score for an insertion of target's AminoAcid n.
	 * 
	 * @param threading
	 *            the threading to be scored
	 * @param n
	 *            the number of the aminoAcid in the target that needs to be
	 *            inserted
	 * @return the score of
	 */
	public double getInsertionScore(Threading threading, int n) {

		// PDBEntry model = modifyModel(threading);
		double result = avgPotential[threading.getSequence().getComp(n) - 65];

		result = -(zeta * 3 * result);
		// result = -1;
		insertmax = insertmax < result ? result : insertmax;
		insertmin = insertmin > result ? result : insertmin;
		insertsum += result;
		insertnum++;
		return result;
	}

	/**
	 * 
	 * @param threading
	 * @param m
	 * @return
	 */
	public double getDeletionScore(Threading threading, int m) {

		// phiP(m)
		// PDBEntry model = modifyModel(threading);
		// double result = pcp.getAminoScores(model)[m];
		double result = avgPotential[threading.getStructure().getAminoAcid(m)
				.getName().getNumber()];

		result = -(zeta * 3 * result);
		// result = - 1;
		deletemax = deletemax < result ? result : deletemax;
		deletemin = deletemin > result ? result : deletemin;
		deletesum += result;
		deletenum++;
		return result;
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

}

/******************************************************************************
 * "A question that sometimes drives me hazy: * Am I or are the others crazy?" *
 * - Albert Einstein (1879 - 1955) *
 ******************************************************************************/
