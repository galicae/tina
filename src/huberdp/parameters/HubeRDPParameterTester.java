/******************************************************************************
 * validation.huberdp.HubeRDPValidator.java                                   *
 *                                                                            *
 * This file is best read at line width 80 and tab width 4.                   *
 *                                                                   huberste *
 ******************************************************************************/
package huberdp.parameters;

import java.text.DecimalFormat;
import java.util.Locale;

import bioinfo.Sequence;
import bioinfo.alignment.SequenceAlignment;
import bioinfo.alignment.matrices.QuasarMatrix;
import bioinfo.energy.potential.SipplContactPotential;
import bioinfo.energy.potential.hydrophobicity.HydrophobicityMatrix;
import bioinfo.proteins.CCPMatrix;
import bioinfo.proteins.PDBEntry;
import bioinfo.proteins.PDBFileReader;
import bioinfo.superpos.TMMain;
import bioinfo.superpos.Transformation;
import huberdp.HubeRDP;
import huberdp.RDPPriorityQueue;
import huberdp.RDPProblem;
import huberdp.RDPSolutionTree;
import huberdp.oracles.*;
import huberdp.scoring.*;
import util.CommandLineParser;

/**
 * HubeRDPParameterOptimzer changes parameters for parameter optimization. use
 * with script & on grid!
 * 
 * @author huberste
 * @lastchange 2013-02-22
 */
public class HubeRDPParameterTester {

	private static boolean test = false;

	private static final String TARGETSTRING = "1dp7P00";
	private static final String TEMPLATESTRING = "1j2xA00";
	private static final String PDBPATH = "/home/h/huberste/gobi/data/pdb/";
	// private static final String GAMMASTRING = "1.0";
	// private static final String DELTASTRING = "0.1";
	// private static final String EPSILONSTRING = "2.0";
	private static final String ZETASTRING = "4.0";
	private static final String GAPSTRING = "3.0";
	private static final String HYDROSTRING = "/home/h/huberste/gobi/data/hydrophobicityMatrices/hydro_1024";
	private static final String CCPSTRING = "/home/h/huberste/gobi/data/CCP/ccp";
	private static final String VOROSTRING = "/home/h/huberste/gobi/tina/tools/voro++_ubuntuquantal";
	private static final String VPOTSTRING = "/home/h/huberste/gobi/data/vpot/s3d3_hob97_25_ED6_SD6_cb_all.md15.hssp95.vpot";
	private static final String DSSPSTRING = "/home/h/huberste/gobi/data/dssp/";
	private static final String TEMPDIR = "/tmp/";

	private static String USAGE = "";

	/**
	 * Validates the HubeRDP.
	 * 
	 * @param args
	 * @throws Exception
	 */
	public static void main(String[] args) {

		// allocate memory
		String targetString = null, templateString = null, pdbpath = null, hydroFileString = null, ccpFileString = null, voroFileString = null, vpotFileString = null, dsspPathString = null, tempdir = null;
		double /* gamma = 0.0, delta = 0.0, epsilon = 0.0, */zeta = 0.0, gap = 0.0;
		if (test) {
			targetString = TARGETSTRING;
			templateString = TEMPLATESTRING;
			pdbpath = PDBPATH;
			hydroFileString = HYDROSTRING;
			ccpFileString = CCPSTRING;
			voroFileString = VOROSTRING;
			vpotFileString = VPOTSTRING;
			dsspPathString = DSSPSTRING;
			// gamma = Double.parseDouble(GAMMASTRING);
			// delta = Double.parseDouble(DELTASTRING);
			// epsilon = Double.parseDouble(EPSILONSTRING);
			zeta = Double.parseDouble(ZETASTRING);
			gap = Double.parseDouble(GAPSTRING);
			tempdir = TEMPDIR;
		} else {
			CommandLineParser clp = new CommandLineParser(args);
			if ((targetString = clp.getStringArg("--target")) == null) {
				System.out.println(USAGE);
				System.out.println("No --target was given");
				System.exit(0);
			}
			if ((templateString = clp.getStringArg("--template")) == null) {
				System.out.println(USAGE);
				System.out.println("No --template was given");
				System.exit(0);
			}
			if ((pdbpath = clp.getStringArg("--pdb")) == null) {
				System.out.println(USAGE);
				System.out.println("No --pdb was given");
				System.exit(0);
			}
			if ((hydroFileString = clp.getStringArg("--hydro")) == null) {
				System.out.println(USAGE);
				System.out.println("No --hydro was given");
				System.exit(0);
			}
			if ((ccpFileString = clp.getStringArg("--ccp")) == null) {
				System.out.println(USAGE);
				System.out.println("No --ccp was given");
				System.exit(0);
			}
			if ((voroFileString = clp.getStringArg("--voro")) == null) {
				System.out.println(USAGE);
				System.out.println("No --voro was given");
				System.exit(0);
			}
			if ((vpotFileString = clp.getStringArg("--vpot")) == null) {
				System.out.println(USAGE);
				System.out.println("No --vpot was given");
				System.exit(0);
			}
			if ((dsspPathString = clp.getStringArg("--dssp")) == null) {
				System.out.println(USAGE);
				System.out.println("No --dssp was given");
				System.exit(0);
			}
			// if ((gamma = clp.getDoubleArg("--gamma")) == 0.0) {
			// System.out.println(usage);
			// System.out.println("No --gamma was given");
			// System.exit(0);
			// }
			// if ((delta = clp.getDoubleArg("--delta")) == 0.0) {
			// System.out.println(usage);
			// System.out.println("No --delta was given");
			// System.exit(0);
			// }
			// if ((epsilon = clp.getDoubleArg("--epsilon")) == 0.0) {
			// System.out.println(usage);
			// System.out.println("No --epsilon was given");
			// System.exit(0);
			// }
			if ((zeta = clp.getDoubleArg("--zeta")) == 0.0) {
				System.out.println(USAGE);
				System.out.println("no --zeta was given");
				System.exit(1);
			}
			if ((gap = clp.getDoubleArg("--gap")) == 0.0) {
				System.out.println(USAGE);
				System.out.println("no --gap was given");
				System.exit(1);
			}
			if ((tempdir = clp.getStringArg("--temp")) == null) {
				System.out.println(USAGE);
				System.out.println("no --temp was given");
				System.exit(1);
			}
			clp = null; // clp -> GC
		}

		// initialize important stuff
		Locale.setDefault(Locale.US);
		DecimalFormat df = new DecimalFormat("0.0000");
		// initialize TM stuff
		TMMain tmmain = new TMMain();
		// initialize PDBFiles
		PDBEntry templateStructure = new PDBFileReader(pdbpath)
				.readPDBFromFile(pdbpath + templateString + ".pdb");
		PDBEntry targetStructure = new PDBFileReader(pdbpath)
				.readPDBFromFile(pdbpath + targetString + ".pdb");
		Sequence targetSequence = targetStructure.getSequence();

		SipplContactPotential sippl = new SipplContactPotential();
		sippl.readFromVPOTFile(vpotFileString);

		HydrophobicityMatrix hydroM = new HydrophobicityMatrix(hydroFileString);
		CCPMatrix ccpM = new CCPMatrix(ccpFileString);

		for (double gammavar = 1; gammavar <= 10; gammavar += 1) {
			for (double deltavar = 1; deltavar <= 10; deltavar += 1) {
				for (double epsilonvar = 1; epsilonvar <= 10; epsilonvar += 1) {

					// align sequences with HubeRDP
					// initialize HubeRDP stuff
					HubeRDP rdp = new HubeRDP();
					RDPScoring scoring = new RDPScoring(gammavar, deltavar,
							epsilonvar, zeta, gap, QuasarMatrix.DAYHOFF_MATRIX,
							hydroM, ccpM, dsspPathString, sippl,
							templateStructure, voroFileString,
							RDPScoring.GRID_EXTEND, RDPScoring.GRID_DENSITY,
							RDPScoring.GRID_CLASH, RDPScoring.MIN_CONTACT,
							tempdir);
					rdp.setScoring(scoring);
					rdp.addOracle(new RDPOracle(scoring));

					RDPProblem root = new RDPProblem(templateStructure,
							targetSequence);
					RDPSolutionTree t = new RDPSolutionTree(root);
					RDPPriorityQueue pq = new RDPPriorityQueue(t.getRoot());
					// run HubeRDP
					rdp.rdp(t, pq);

					// get HubeRDP's (first) alignment
					SequenceAlignment rdpAlignment = t.getRoot().getTA().get(0)
							.getThreading().asSequenceAlignment();

					Transformation rdptmtr = null;
					try {
						rdptmtr = tmmain.calculateTransformation(rdpAlignment,
								templateStructure, targetStructure);
						System.out.println(df.format(gammavar)
								+ " "
								+ df.format(deltavar)
								+ " "
								+ df.format(epsilonvar)
								+ " "
								+ df.format(zeta)
								+ " "
								+ df.format(gap)
								+ " "
								+ df.format(t.getRoot().getTA().getFirst()
										.getThreading().getScore()) + " "
								+ df.format(rdptmtr.getRmsd()) + " "
								+ df.format(rdptmtr.getGdt()) + " "
								+ df.format(rdptmtr.getTmscore()));
					} catch (Exception e) {
						System.err.println("error calculating TM Stuff: "
								+ e.getLocalizedMessage());
						e.printStackTrace();
					}
				}
			}
		}
	}
}

/******************************************************************************
 * "A question that sometimes drives me hazy: * Am I or are the others crazy?" *
 * - Albert Einstein (1879 - 1955) *
 ******************************************************************************/
