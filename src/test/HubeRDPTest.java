/******************************************************************************
 * test.HubeRDPTest.java                                                      *
 *                                                                            *
 * Contains a main method for calling and testing HubeRDP.                    *
 *                                                                            *
 * This file is best read at line width 80 and tab width 4.                   *
 *                                                                   huberste *
 ******************************************************************************/
package test;

import huberdp.HubeRDP;
import huberdp.RDPPriorityQueue;
import huberdp.RDPProblem;
import huberdp.RDPSolutionTree;
import huberdp.Scoring;
import huberdp.oracles.RDPOracle;
import huberdp.scoring.MagicScoring;
import huberdp.scoring.RDPScoring;
import huberdp.scoring.SimpleScoring;

import java.text.DecimalFormat;
import java.util.Locale;

import bioinfo.alignment.SequenceAlignment;
import bioinfo.alignment.gotoh.Gotoh;
import bioinfo.alignment.gotoh.LocalSequenceGotoh;
import bioinfo.alignment.matrices.QuasarMatrix;
import bioinfo.energy.potential.SipplContactPotential;
import bioinfo.energy.potential.hydrophobicity.HydrophobicityMatrix;
import bioinfo.proteins.CCPMatrix;
import bioinfo.proteins.PDBEntry;
import bioinfo.proteins.PDBFileReader;
import bioinfo.superpos.TMMain;
import bioinfo.superpos.Transformation;

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
	 * main function for calling HubeRDP
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
		Scoring simpleScoring = new SimpleScoring();
		Scoring magicScoring = new MagicScoring(PDBSTRING, TMALIGN, TEMPDIR);

		SipplContactPotential sippl = new SipplContactPotential();
		sippl.readFromVPOTFile(VPOTSTRING);
		Scoring rdpScoring = new RDPScoring(GAMMA, DELTA, EPSILON, ZETA, GAP,
				QuasarMatrix.DAYHOFF_MATRIX, new HydrophobicityMatrix(
						HYDROSTRING), new CCPMatrix(CCPSTRING), DSSPSTRING,
				sippl, templateStructure, VOROSTRING, RDPScoring.GRID_EXTEND,
				RDPScoring.GRID_DENSITY, RDPScoring.GRID_CLASH,
				RDPScoring.MIN_CONTACT, TEMPDIR);

		Scoring scoring = rdpScoring;

		rdp.setScoring(scoring);

		// add oracles
		rdp.addOracle(new RDPOracle(scoring));

		// execute rdp algorithm
		rdp.rdp(t, pq);
		// Solutions are now in t.getRoot();

		// get HubeRDP's (first) alignment
		SequenceAlignment rdpAlignment = t.getRoot().getTA().get(0)
				.getThreading().asSequenceAlignment();

		// initialize TM stuff
		TMMain tmmain = new TMMain();

		Transformation rdptmtr = tmmain.calculateTransformation(rdpAlignment,
				templateStructure, targetStructure);

		// initialize output stuff
		Locale.setDefault(Locale.US);
		DecimalFormat df = new DecimalFormat("0.0000");

		System.out.println("RDP:");
		System.out.println(rdpAlignment.toStringVerbose());
		System.out.println("RDP-Scr\tRMSD\tGDT\tTM-Scr");
		System.out.println(df.format(rdpAlignment.getScore()) + "\t"
				+ df.format(rdptmtr.getRmsd()) + "\t"
				+ df.format(rdptmtr.getGdt()) + "\t"
				+ df.format(rdptmtr.getTmscore()));

		Gotoh gotoh = new LocalSequenceGotoh(-10.0, -2.0, QuasarMatrix.DAYHOFF_MATRIX);
		SequenceAlignment gotAlignment = (SequenceAlignment) gotoh.align(targetStructure.getSequence(), templateStructure.getSequence());
		Transformation gottmtr = tmmain.calculateTransformation(gotAlignment,
				templateStructure, targetStructure);
		
		System.out.println("Gotoh:");
		System.out.println(gotAlignment.toStringVerbose());
		System.out.println("gth-Scr\tRMSD\tGDT\tTM-Scr");
		System.out.println(df.format(gotAlignment.getScore()) + "\t"
				+ df.format(gottmtr.getRmsd()) + "\t"
				+ df.format(gottmtr.getGdt()) + "\t"
				+ df.format(gottmtr.getTmscore()));
		
		
		// PREDICT STRUCTURE
		// use CoordMapper on HubeRDP Alignment to construct Structure from
		// Alignment
		// PDBEntry rdpStructure = SimpleCoordMapper.map(templateStructure,
		// rdpAlignment);
		// rdpStructure.writeToPDBFile("/home/h/huberste/gobi/out/"+templateStructure.getID()+targetStructure.getID()+".pdb");

	}

}

/******************************************************************************
 * "A question that sometimes drives me hazy: * Am I or are the others crazy?" *
 * - Albert Einstein (1879 - 1955) *
 ******************************************************************************/
