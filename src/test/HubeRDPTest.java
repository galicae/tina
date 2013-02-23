/******************************************************************************
 * test.HubeRDPTest.java                                                      *
 * huberdptest is mainly a test routine for calling hubeRDP.                  *
 *                                                                            *
 * This file is best read at line width 80 and tab width 4.                   *
 *                                                                   huberste *
 ******************************************************************************/
package test;

import java.text.DecimalFormat;
import java.util.Locale;

import bioinfo.alignment.SequenceAlignment;
import bioinfo.alignment.matrices.QuasarMatrix;
import bioinfo.energy.potential.SipplContactPotential;
import bioinfo.energy.potential.hydrophobicity.HydrophobicityMatrix;
import bioinfo.proteins.CCPMatrix;
import bioinfo.proteins.PDBEntry;
import bioinfo.proteins.PDBFileReader;
import bioinfo.proteins.structure.SimpleCoordMapper;
import bioinfo.superpos.TMMain;
import bioinfo.superpos.Transformation;
import huberdp.*;
import huberdp.oracles.*;
import huberdp.scoring.*;

/**
 * 
 * @author huberste
 * @lastchange 2013-02-22
 */
public class HubeRDPTest {

	private static final String TEMPLATESTRING = "/home/h/huberste/gobi/data/pdb/1j2xA00.pdb";
	private static final String TARGETSTRING = "/home/h/huberste/gobi/data/pdb/1dp7P00.pdb";
	private static final String PDBSTRING = "/home/h/huberste/gobi/data/pdb/";
	private static final String TMALIGN = "/home/h/huberste/gobi/tina/tools/TMalign";
	private static final String HYDROSTRING = "/home/h/huberste/gobi/data/hydrophobicityMatrices/hydro_1024";
	private static final String CCPSTRING = "/home/h/huberste/gobi/data/CCP/ccp";
	private static final String VOROSTRING = "/home/h/huberste/gobi/tina/tools/voro++_ubuntuquantal";
	private static final String VPOTSTRING = "/home/h/huberste/gobi/data/vpot/s3d3_hob97_25_ED6_SD6_cb_all.md15.hssp95.vpot";
	private static final String DSSPSTRING = "/home/h/huberste/gobi/data/dssp/";
	private static final String TEMPDIR = "/tmp/";

	private static final double GAMMA = 1.0, DELTA = 0.1, EPSILON = 2.0,
			ZETA = 4.0, GAP = 14.0;

	/**
	 * main function
	 * 
	 * @param args
	 *            no args needed
	 * @throws Exception
	 */
	public static void main(String[] args) throws Exception {

		// set test data
		PDBEntry templateStructure = null;
		PDBEntry targetStructure = null;

		// read teplate pdb file
		PDBFileReader fr = new PDBFileReader();
		templateStructure = fr.readPDBFromFile(TEMPLATESTRING);
		targetStructure = fr.readPDBFromFile(TARGETSTRING);
		// nullify fr for the Garbage Collector
		fr = null;

		// construct rdp tree
		RDPProblem root = new RDPProblem(templateStructure,
				targetStructure.getSequence());
		RDPSolutionTree t = new RDPSolutionTree(root);

		// construct priority queue
		RDPPriorityQueue pq = new RDPPriorityQueue(t.getRoot());

		// construct RDP
		HubeRDP rdp = new HubeRDP();

		// set scoring
		// rdp.setScoring(new ManualScoring());
		// rdp.setScoring(new SimpleScoring());
		SipplContactPotential sippl = new SipplContactPotential();
		sippl.readFromVPOTFile(VPOTSTRING);
		RDPScoring scoring = new RDPScoring(GAMMA, DELTA, EPSILON, ZETA, GAP,
				QuasarMatrix.DAYHOFF_MATRIX, new HydrophobicityMatrix(
						HYDROSTRING), new CCPMatrix(CCPSTRING), DSSPSTRING,
				sippl, templateStructure, VOROSTRING, RDPScoring.GRID_EXTEND,
				RDPScoring.GRID_DENSITY, RDPScoring.GRID_CLASH,
				RDPScoring.MIN_CONTACT, TEMPDIR);
		rdp.setScoring(scoring);
//		rdp.setScoring(new SimpleScoring());
//		MagicScoring magic = new MagicScoring(PDBSTRING, TMALIGN);
//		rdp.setScoring(magic);

		// add oracles
		// rdp.addOracle(new ManualOracle());
		// rdp.addOracle(new TinyOracle());
		rdp.addOracle(new RDPOracle(scoring));
//		rdp.addOracle(new TinyOracle());

		// execute rdp algorithm
		rdp.rdp(t, pq);
		// Solutions are now in t.getRoot();

		// get HubeRDP's (first) alignment
		SequenceAlignment rdpAlignment = t.getRoot().getTA().get(0)
				.getThreading().asSequenceAlignment();

		// SequenceAlignment rdpAlignment = new SequenceAlignment(
		// templateStructure.getSequence(),
		// targetStructure.getSequence(),
		// "GPLDVQVT--EDAVRRYLTR-KPMTTKDLLKKF-Q-TKKTGLSSEQTVNVLAQILKR-LNPERKMIN----DKMHFSLKE----",
		// "------TVQWLLDNYE-TAEGVSLPRSTLYNHYLLHSQEQKLEPVNAASFGKLIRSVFMGLRTRRLGTRGNSKYHYYGL-RIKA",
		// 0.0);

		System.out.println(rdpAlignment.toStringVerbose());

		// initialize TM stuff
		TMMain tmmain = new TMMain();

		Transformation rdptmtr = tmmain.calculateTransformation(rdpAlignment,
				templateStructure, targetStructure);

		// initialize output stuff
		Locale.setDefault(Locale.US);
		DecimalFormat df = new DecimalFormat("0.0000");
		DecimalFormat dfshort = new DecimalFormat("0.00");

		System.out
				.println("gamma\tdelta\tepsilon\tzeta\tgap\tRDP-Scr\tRMSD\tGDT\tTM-Scr");
		System.out.println(dfshort.format(GAMMA)
				+ "\t"
				+ dfshort.format(DELTA)
				+ "\t"
				+ dfshort.format(EPSILON)
				+ "\t"
				+ dfshort.format(ZETA)
				+ "\t"
				+ dfshort.format(GAP)
				+ "\t"
				+ df.format(t.getRoot().getTA().getFirst().getThreading()
						.getScore()) + "\t" + df.format(rdptmtr.getRmsd())
				+ "\t" + df.format(rdptmtr.getGdt()) + "\t"
				+ df.format(rdptmtr.getTmscore()));

		// PREDICT STRUCTURE
		// use CoordMapper on HubeRDP Alignment to construct Structure from
		// Alignment
//		 PDBEntry rdpStructure = SimpleCoordMapper.map(templateStructure,
//		 rdpAlignment);
//		 rdpStructure.writeToPDBFile("/home/h/huberste/gobi/out/"+templateStructure.getID()+targetStructure.getID()+".pdb");

	}

}

/******************************************************************************
 * "A question that sometimes drives me hazy: * Am I or are the others crazy?" *
 * - Albert Einstein (1879 - 1955) *
 ******************************************************************************/
