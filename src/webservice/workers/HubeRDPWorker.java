/**
 * 
 */
package webservice.workers;

import huberdp.HubeRDP;
import huberdp.RDPPriorityQueue;
import huberdp.RDPProblem;
import huberdp.RDPSolutionTree;
import huberdp.Scoring;
import huberdp.oracles.RDPOracle;
import huberdp.scoring.MagicScoring;
import huberdp.scoring.RDPScoring;
import huberdp.scoring.SimpleScoring;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;

import bioinfo.alignment.Threading;
import bioinfo.alignment.matrices.QuasarMatrix;
import bioinfo.energy.potential.SipplContactPotential;
import bioinfo.energy.potential.hydrophobicity.HydrophobicityMatrix;
import bioinfo.proteins.CCPMatrix;
import bioinfo.proteins.PDBEntry;
import bioinfo.proteins.PDBFileReader;

/**
 * @author huberste
 * @lastchange 2013-02-26
 */
public class HubeRDPWorker extends Worker {

	private static final String PDBSTRING = "/home/h/huberste/gobi/data/pdb/";
	private static final String HYDROSTRING = "/home/h/huberste/gobi/data/hydrophobicityMatrices/hydro_1024";
	private static final String CCPSTRING = "/home/h/huberste/gobi/data/CCP/ccp";
	private static final String VOROSTRING = "/home/h/huberste/gobi/tina/tools/voro++_ubuntuquantal";
	private static final String VPOTSTRING = "/home/h/huberste/gobi/data/vpot/vpot";
	private static final String DSSPSTRING = "/home/h/huberste/gobi/data/dssp/";
	private static final String TEMPDIR = "/tmp/";

	private static final double GAMMA = 1.0, DELTA = 0.1, EPSILON = 2.0,
			ZETA = 4.0, GAP = 14.0;

	private String template;
	private String target;

	private Threading result;

	public HubeRDPWorker(String jobFile) {
		super(jobFile);
	}

	@Override
	public void work() {
		// read WORKING job file
		readFile();

		// set test data
		PDBEntry templateStructure = null;
		PDBEntry targetStructure = null;

		// read teplate pdb file
		PDBFileReader fr = new PDBFileReader();
		templateStructure = fr.readPDBFromFile(PDBSTRING + template + ".pdb");
		targetStructure = fr.readPDBFromFile(PDBSTRING + target + ".pdb");
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
		SipplContactPotential sippl = new SipplContactPotential();
		sippl.readFromVPOTFile(VPOTSTRING);
		Scoring scoring = new RDPScoring(GAMMA, DELTA, EPSILON, ZETA, GAP,
				QuasarMatrix.DAYHOFF_MATRIX, new HydrophobicityMatrix(
						HYDROSTRING), new CCPMatrix(CCPSTRING), DSSPSTRING,
				sippl, templateStructure, VOROSTRING, RDPScoring.GRID_EXTEND,
				RDPScoring.GRID_DENSITY, RDPScoring.GRID_CLASH,
				RDPScoring.MIN_CONTACT, TEMPDIR);

		rdp.setScoring(scoring);

		// add oracles
		rdp.addOracle(new RDPOracle(scoring));

		// execute rdp algorithm
		rdp.rdp(t, pq);
		// Solutions are now in t.getRoot();

		// get HubeRDP's (first) alignment
		result = t.getRoot().getTA().get(0).getThreading();

		// Write DONE job file
		writeResult();
	}

	@Override
	protected void readFile() {
		BufferedReader from = null;

		String line = null;
		try {
			from = new BufferedReader(new FileReader(JOB_FILE));
			while ((line = from.readLine()) != null) {
				if (line.startsWith("TEMPLATE_ID=")) {
					String[] temp = line.substring(13).split("=");
					template = temp[1];
				} else if (line.startsWith("TARGET_ID=")) {
					String[] temp = line.substring(13).split("=");
					target = temp[1];
				}
			}
		} catch (IOException e) {
			System.err.println("Error while trying to read " + JOB_FILE + ".");
			e.printStackTrace();
		} finally {
			try {
				from.close();
			} catch (IOException e) {
				System.err
						.println("Error while trying close " + JOB_FILE + ".");
				e.printStackTrace();
			}
		}
	}

	@Override
	protected void writeResult() {
		BufferedReader from = null;
		BufferedWriter to = null;

		String line = null;
		try {
			from = new BufferedReader(new FileReader(JOB_FILE));
			to = new BufferedWriter(new FileWriter(DONE_FILE));
			while ((line = from.readLine()) != null) {
				to.write(line + "\n");
			}
			to.write("RESULT=\n");
			to.write(result.toStringVerbose());
		} catch (IOException e) {
			System.err.println("Error while trying to copy " + JOB_FILE
					+ " to " + DONE_FILE + ".");
			e.printStackTrace();
		} finally {
			try {
				if (from != null)
					from.close();
				if (to != null)
					to.close();
			} catch (IOException e) {
				System.err.println("Error while trying close FileStreams");
				e.printStackTrace();
			}
		}
		new File(JOB_FILE).delete();
	}

}
