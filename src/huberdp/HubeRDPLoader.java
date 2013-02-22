/******************************************************************************
 * huberdp.HubeRDPLoader.java                                                 *
 *                                                                            *
 * This file contains the main method to call the HubeRDP as an Threader.     *
 *                                                                            *
 * This file is best read at line width 80 and tab width 4.                   *
 *                                                                   huberste *
 ******************************************************************************/
package huberdp;

import java.io.File;

import static bioinfo.Sequence.getSequenceFromString;

import huberdp.oracles.RDPOracle;
import huberdp.scoring.RDPScoring;
import util.CommandLineParser;
import bioinfo.Sequence;
import bioinfo.alignment.matrices.QuasarMatrix;
import bioinfo.energy.potential.SipplContactPotential;
import bioinfo.energy.potential.hydrophobicity.HydrophobicityMatrix;
import bioinfo.proteins.CCPMatrix;
import bioinfo.proteins.PDBEntry;
import bioinfo.proteins.PDBFileReader;

/**
 * @author huberste
 * @lastchange 2013-02-20
 */
public class HubeRDPLoader {

	private static final boolean test = false;

	private static final String TARGETESTRING = "1dp7P00:TVQWLLDNYETAEGVSLPRSTLYNHYLLHSQEQKLEPVNAASFGKLIRSVFMGLRTRRLGTRGNSKYHYYGLRIK";
	private static final String TEMPLATESTRING = "/home/h/huberste/gobi/data/pdb/1j2xA00.pdb";
	private static final String GAMMASTRING = "1.0";
	private static final String DELTASTRING = "0.1";
	private static final String EPSILONSTRING = "2.0";
	private static final String ZETASTRING = "4.0";
	private static final String GAPSTRING = "3.0";
	private static final String HYDROSTRING = "/home/h/huberste/gobi/data/hydrophobicityMatrices/hydro_1024";
	private static final String CCPSTRING = "/home/h/huberste/gobi/data/CCP/ccp";
	private static final String VOROSTRING = "/home/h/huberste/gobi/tina/tools/voro++_ubuntuquantal";
	private static final String VPOTSTRING = "/home/h/huberste/gobi/data/vpot/s3d3_hob97_25_ED6_SD6_cb_all.md15.hssp95.vpot";
	private static final String DSSPSTRING = "/home/h/huberste/gobi/data/dssp/";
	private static final String TEMPDIR = "/tmp/";

	private static String usage = "HubeRDP Threading Tool\n"
			+ "usage: java -jar HubeRDP.jar \n" + "\t--target <sequence>\n"
			+ "\t--template <templatefile>\n" + "\t--gamma <gamma>\n"
			+ "\t--delta <delta>\n" + "\t--epsilon <epsilon>\n"
			+ "\t--zeta <zeta>\n" + "\t--gap <gap>\n"
			+ "\t--hydro <hydrofile>\n" + "\t--ccp <ccpfile>\n"
			+ "\t--voro <voro++path>\n" + "\t--vpot <vpotfile>\n"
			+ "\t--dssp <dssppath>\n";

	// + "<sequence> is the target Sequence in format \"id:sequence\" and "
	// + "<templatefile> is the (full) path to the PDBFile containing the "
	// + "template's pdbFile. The filename must be in format "
	// + "\"xxxxA00.pdb\" where xxxx is the PDB ID, A is the chain "
	// + "identifier and 00 the number of the model to be used.\n";
	/**
	 * @param args
	 */
	public static void main(String[] args) {
		// allocate memory
		String targetString = null, templateString = null, hydroFileString = null, ccpFileString = null, voroFileString = null, vpotFileString = null, dsspPathString = null, tempdir = null;
		Sequence target = null;
		double gamma = 0.0, delta = 0.0, epsilon = 0.0, zeta = 0.0, gap = 0.0;
		if (test) {
			targetString = TARGETESTRING;
			templateString = TEMPLATESTRING;
			hydroFileString = HYDROSTRING;
			ccpFileString = CCPSTRING;
			voroFileString = VOROSTRING;
			vpotFileString = VPOTSTRING;
			dsspPathString = DSSPSTRING;
			target = getSequenceFromString(targetString);
			gamma = Double.parseDouble(GAMMASTRING);
			delta = Double.parseDouble(DELTASTRING);
			epsilon = Double.parseDouble(EPSILONSTRING);
			zeta = Double.parseDouble(ZETASTRING);
			gap = Double.parseDouble(GAPSTRING);
			tempdir = TEMPDIR;
		} else {
			CommandLineParser clp = new CommandLineParser(args);
			if ((targetString = clp.getStringArg("--target")) == null) {
				System.out.println(usage);
				System.out.println("No --target was given");
				System.exit(0);
			}
			if ((templateString = clp.getStringArg("--template")) == null) {
				System.out.println(usage);
				System.out.println("No --template was given");
				System.exit(0);
			}
			if ((target = getSequenceFromString(targetString)) == null) {
				System.out.println(usage);
				System.out.println("<target> has wrong format");
				System.exit(0);
			}
			if (!(new File(templateString).canRead())) {
				System.out.println(usage);
				System.out.println("can't read template");
				System.exit(0);
			}
			if ((hydroFileString = clp.getStringArg("--hydro")) == null) {
				System.out.println(usage);
				System.out.println("No --hydro was given");
				System.exit(0);
			}
			if ((ccpFileString = clp.getStringArg("--ccp")) == null) {
				System.out.println(usage);
				System.out.println("No --ccp was given");
				System.exit(0);
			}
			if ((voroFileString = clp.getStringArg("--voro")) == null) {
				System.out.println(usage);
				System.out.println("No --voro was given");
				System.exit(0);
			}
			if ((vpotFileString = clp.getStringArg("--vpot")) == null) {
				System.out.println(usage);
				System.out.println("No --vpot was given");
				System.exit(0);
			}
			if ((dsspPathString = clp.getStringArg("--dssp")) == null) {
				System.out.println(usage);
				System.out.println("No --dssp was given");
				System.exit(0);
			}
			if ((gamma = clp.getDoubleArg("--gamma")) == 0.0) {
				System.out.println(usage);
				System.out.println("No --gamma was given");
				System.exit(0);
			}
			if ((delta = clp.getDoubleArg("--delta")) == 0.0) {
				System.out.println(usage);
				System.out.println("No --delta was given");
				System.exit(0);
			}
			if ((epsilon = clp.getDoubleArg("--epsilon")) == 0.0) {
				System.out.println(usage);
				System.out.println("No --epsilon was given");
				System.exit(0);
			}
			if ((zeta = clp.getDoubleArg("--zeta")) == 0.0) {
				System.out.println(usage);
				System.out.println("no --zeta was given");
				System.exit(1);
			}
			if ((gap = clp.getDoubleArg("--gap")) == 0.0) {
				System.out.println(usage);
				System.out.println("no --gap was given");
				System.exit(1);
			}
			if ((tempdir = clp.getStringArg("--temp")) == null) {
				System.out.println(usage);
				System.out.println("no --temp was given");
				System.exit(1);
			}
			clp = null;
		}

		// read template pdb file
		PDBEntry template = new PDBFileReader().readPDBFromFile(templateString);

		// construct rdp tree
		RDPProblem root = new RDPProblem(template, target);
		RDPSolutionTree t = new RDPSolutionTree(root);
		RDPPriorityQueue pq = new RDPPriorityQueue(t.getRoot());
		// construct RDP
		HubeRDP rdp = new HubeRDP();

		// set scoring
		SipplContactPotential sippl = new SipplContactPotential();
		sippl.readFromVPOTFile(vpotFileString);
		RDPScoring scoring = new RDPScoring(gamma, delta, epsilon, zeta, gap,
				QuasarMatrix.DAYHOFF_MATRIX, new HydrophobicityMatrix(
						hydroFileString), new CCPMatrix(ccpFileString),
				dsspPathString, sippl, template, voroFileString,
				RDPScoring.GRID_EXTEND, RDPScoring.GRID_DENSITY,
				RDPScoring.GRID_CLASH, RDPScoring.MIN_CONTACT, tempdir);

		rdp.setScoring(scoring);

		// add oracle(s)
		rdp.addOracle(new RDPOracle(scoring));
		// execute rdp algorithm
		rdp.rdp(t, pq);
		// Solution is now in t.getRoot();
		System.out.println(t.getRoot().getTA().toString());

	}

}

/******************************************************************************
 * "A question that sometimes drives me hazy: * Am I or are the others crazy?" *
 * - Albert Einstein (1879 - 1955) *
 ******************************************************************************/
