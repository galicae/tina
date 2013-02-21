/******************************************************************************
 * test.HubeRDPTest.java                                                      *
 * huberdptest is mainly a test routine for calling hubeRDP.                  *
 *                                                                            *
 * This file is best read at line width 80 and tab width 4.                   *
 *                                                                   huberste *
 ******************************************************************************/
package test;

import bioinfo.Sequence;
import bioinfo.alignment.matrices.QuasarMatrix;
import bioinfo.energy.potential.SipplContactPotential;
import bioinfo.energy.potential.hydrophobicity.HydrophobicityMatrix;
import bioinfo.proteins.CCPMatrix;
import bioinfo.proteins.PDBEntry;
import bioinfo.proteins.PDBFileReader;
import huberdp.*;
import huberdp.oracles.*;
import huberdp.scoring.*;

/**
 * 
 * @author huberste
 * @lastchange 2013-02-14
 */
public class HubeRDPTest {

	private static final String TARGETESTRING = "1dp7P00:TVQWLLDNYETAEGVSLPRSTLYNHYLLHSQEQKLEPVNAASFGKLIRSVFMGLRTRRLGTRGNSKYHYYGLRIK";
	private static final Sequence TARGET = new Sequence("1dp7P00","TVQWLLDNYETAEGVSLPRSTLYNHYLLHSQEQKLEPVNAASFGKLIRSVFMGLRTRRLGTRGNSKYHYYGLRIK");
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
	
	/**
	 * main function
	 * 
	 * @param args
	 *            no args needed
	 */
	public static void main(String[] args) {

		// set test data
		PDBEntry templateStructure = null;

		// read teplate pdb file
		System.out.print("Reading Template structure file...");
		PDBFileReader fr = new PDBFileReader();
		templateStructure = fr
				.readPDBFromFile(TEMPLATESTRING);
		// nullify fr for the Garbage Collector
		fr = null;
		System.out.println(" done!");

		// construct rdp tree
		System.out.print("Constructing RDP Tree structure...");

		RDPProblem root = new RDPProblem(templateStructure, TARGET);
		RDPSolutionTree t = new RDPSolutionTree(root);
		System.out.println(" done!");

		// construct priority queue
		System.out.print("constructing Priority Queue...");
		RDPPriorityQueue pq = new RDPPriorityQueue(t.getRoot());
		System.out.println(" done!");

		// construct RDP
		HubeRDP rdp = new HubeRDP();

		// set scoring
		// rdp.setScoring(new ManualScoring());
		// rdp.setScoring(new SimpleScoring());
		SipplContactPotential sippl = new SipplContactPotential();
		sippl.readFromVPOTFile(VPOTSTRING);
		RDPScoring scoring = new RDPScoring(RDPScoring.GAMMA, RDPScoring.DELTA,
				RDPScoring.EPSILON, RDPScoring.ZETA, RDPScoring.GAP,
				QuasarMatrix.DAYHOFF_MATRIX, new HydrophobicityMatrix(
						HYDROSTRING), new CCPMatrix(CCPSTRING),
				"/home/h/huberste/gobi/data/dssp/", sippl, templateStructure,
				RDPScoring.VOROPATH, RDPScoring.GRID_EXTEND,
				RDPScoring.GRID_DENSITY, RDPScoring.GRID_CLASH,
				RDPScoring.MIN_CONTACT);
		rdp.setScoring(scoring);

		// add oracles
		// rdp.addOracle(new ManualOracle());
		// rdp.addOracle(new TinyOracle());
		rdp.addOracle(new RDPOracle(scoring));

		// execute rdp algorithm
		System.out.println("HubeRDP will now be executed!");
		rdp.rdp(t, pq);
		System.out.println("HubeRDP was successfully executed!");
		// Solution is now in t.getRoot();

		System.out.println(t.getRoot().getTA().toString());

	}
}

/******************************************************************************
 * "A question that sometimes drives me hazy: * Am I or are the others crazy?" *
 * - Albert Einstein (1879 - 1955) *
 ******************************************************************************/
