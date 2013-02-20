/******************************************************************************
 * test.HubeRDPTest.java                                                      *
 * huberdptest is mainly a test routine for calling hubeRDP.                  *
 *                                                                            *
 * This file is best read at line width 80 and tab width 4.                   *
 *                                                                   huberste *
 ******************************************************************************/
package test;

import bioinfo.Sequence;
import bioinfo.alignment.SequenceAlignment;
import bioinfo.alignment.Threading;
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

	private static final String hydromatrixFile = "/home/h/huberste/gobi/data/hydrophobicityMatrices/hydro_1024";
	private static final String ccpfilename = "/home/h/huberste/gobi/data/ccp/ccp";
	private static final String dsspFile = "/home/h/huberste/gobi/data/dssp/";
	private static final String vpotfile = "/home/h/huberste/gobi/data/pair/s3d3_hob97_25_ED6_SD6_cb_all.md15.hssp95.vpot";

	// private static final Sequence template = new
	// Sequence("test001","GGGGCA");
	// private static final Sequence target = new
	// Sequence("test002","TTTGGGGA");

	/**
	 * main function
	 * 
	 * @param args
	 *            no args needed
	 */
	public static void main(String[] args) {

		// set test data
		// set test data
		Sequence template = new Sequence("1j2xA00",
				"GPLDVQVTEDAVRRYLTRKPMTTKDLLKKFQTKKTGLSSEQTVNVLAQILKRLNPERKMINDKMHFSLK");
		Sequence target = new Sequence("1dp7P00",
				"TVQWLLDNYETAEGVSLPRSTLYNHYLLHSQEQKLEPVNAASFGKLIRSVFMGLRTRRLGTRGNSKYHYYGLRIK");
		PDBEntry templateStructure = null;
		PDBEntry targetStructure = null;
		SequenceAlignment ali = null;

		// read teplate pdb file
		System.out.print("Reading Template structure file...");
		PDBFileReader fr = new PDBFileReader();
		templateStructure = fr
				.readPDBFromFile("/home/h/huberste/gobi/webserver/pdb/1J2XA00.pdb");
		// nullify fr for the Garbage Collector
		fr = null;
		template = templateStructure.getSequence();
		System.out.println(" done!");

		// construct rdp tree
		System.out.print("Constructing RDP Tree structure...");
		int[][] rows = new int[2][templateStructure.length() + target.length()];
		for(int i = 0; i < template.length(); i++) {
			rows[0][i] = i;
			rows[1][i] = -1;
		}
		for(int i = 0; i < target.length(); i++) {
			rows[0][templateStructure.length()+i] = -1;
			rows[1][templateStructure.length()+i] = i;
		}
		
		RDPProblem root = new RDPProblem(new Threading(templateStructure, target, rows, 0), 0, rows[0].length-1);
		RDPSolutionTree t = new RDPSolutionTree(root);
		System.out.println(" done!");

		// construct priority queue
		System.out.print("constructing Priority Queue...");
		RDPPriorityQueue pq = new RDPPriorityQueue(t.getRoot());
		System.out.println(" done!");

		// construct RDP
		HubeRDP rdp = new HubeRDP();

		// add oracles
//		rdp.addOracle(new ManualOracle());
//		rdp.addOracle(new TinyOracle());
		rdp.addOracle(new RDPOracle());

		// set scoring
//		rdp.setScoring(new ManualScoring());
//		rdp.setScoring(new SimpleScoring());
		SipplContactPotential sippl = new SipplContactPotential();
		sippl.readFromVPOTFile(vpotfile);
		rdp.setScoring(new RDPScoring(RDPScoring.GAMMA, RDPScoring.DELTA,
				RDPScoring.EPSILON, RDPScoring.ZETA,
				QuasarMatrix.DAYHOFF_MATRIX, new HydrophobicityMatrix(
						hydromatrixFile), new CCPMatrix(ccpfilename), sippl,
				templateStructure, RDPScoring.VOROPATH, RDPScoring.GRID_EXTEND,
				RDPScoring.GRID_DENSITY, RDPScoring.GRID_CLASH,
				RDPScoring.MIN_CONTACT));

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
