/******************************************************************************
 * validation.huberdp.HubeRDPValidator.java                                   *
 *                                                                            *
 * This file is best read at line width 80 and tab width 4.                   *
 *                                                                   huberste *
 ******************************************************************************/
package validation.huberdp;

import java.text.DecimalFormat;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.Locale;

import bioinfo.Sequence;
import bioinfo.alignment.SequenceAlignment;
import bioinfo.alignment.Threading;
import bioinfo.alignment.gotoh.LocalSequenceGotoh;
import bioinfo.alignment.matrices.QuasarMatrix;
import bioinfo.energy.potential.SipplContactPotential;
import bioinfo.energy.potential.hydrophobicity.HydrophobicityMatrix;
import bioinfo.pdb.PDBFile;
import bioinfo.proteins.CCPMatrix;
import bioinfo.proteins.PDBEntry;
import bioinfo.proteins.PDBFileReader;
import bioinfo.proteins.structure.SimpleCoordMapper;
import bioinfo.superpos.TMMain;
import bioinfo.superpos.Transformation;
import files.PairFile;
import files.SeqlibFile;
import huberdp.HubeRDP;
import huberdp.RDPPriorityQueue;
import huberdp.RDPProblem;
import huberdp.RDPSolutionTree;
import huberdp.oracles.*;
import huberdp.scoring.*;
import util.CommandLineParser;

/**
 * HubeRDPValidator validates HubeRDP as an Sequence Alignment Tool against
 * Gotoh.
 * @author huberste
 * @lastchange 2013-02-14
 */
public class HubeRDPValidator {

	private static final String hydromatrixFile = "/home/h/huberste/gobi/data/hydrophobicityMatrices/hydro_1024";
	private static final String ccpfilename = "/home/h/huberste/gobi/data/ccp/ccp";
	private static final String vpotfile = "/home/h/huberste/gobi/data/pair/s3d3_hob97_25_ED6_SD6_cb_all.md15.hssp95.vpot";

	private static String usage =
	"Alignment aller Paare in pairfile:\n"+
	"java HubeRDPValidator --seqlib <seqlibfile> --pairs <pairfile> "+
		"\t--pdbpath <pdbfilepath>\n\n"+
	"<seqlibfile> enthaelt Zeilen der Form \"id:sequence\"\n"+
	"<pairfile> enthaelt Zeilen der Form \"idone idtwo\"\n\n"+
	"<pdbfilepath> is the path where the pdbFiles are in."+
/*
	"Optionale Parameter sind:\n"+
	"\t-m\t\tmatrixname (Standard dayhoff)\n"+
	"\t-go\t\tgapopen (Standard -12)\n"+
	"\t-ge\t\tgapextend (Standard -1)\n"+
	"\t-mode\t\teines aus local|global|freeshift (Standard freeshift)\n"+
	"\t-printali\tgibt auch jedes Alignment aus\n"+
	"\t-printmatrices\teines aus txt|html, gibt auch die Gotoh-Matrizen aus,\n"+
	"\t\t\tentweder als tab separiert oder html\n"+
	"\t-check\t\tueberprueft die berechneten Scores anhand des Alignments\n\n"+
*/
	"Beispiel-Aufruf:\n"+
	"java HubeRDPValidator --seqlib domains.seqlib --pairs cathscop.inpairs \n"+
		"\t--pdbpath";
	
	/**
	 * Validates the HubeRDP.
	 * @param args
	 * @throws Exception 
	 */
	public static void main(String[] args) throws Exception {
		
		// get command line arguments
		CommandLineParser clp = new CommandLineParser(args);
		String seqlib = clp.getStringArg("--seqlib");
		String pairs = clp.getStringArg("--pairs");
		String pdbpath = clp.getStringArg("--pdbpath");
		
		if (pairs == null || seqlib == null || pdbpath == null) {
			System.out.println(usage);
			System.exit(1);
		}
		
		// initialize important stuff	
		Locale.setDefault(Locale.US);
		DecimalFormat df = new DecimalFormat("0.0000");
		
		// load joblist
		PairFile pairfile = new PairFile(pairs);
		LinkedList<String[]> joblist = pairfile.getJoblist();
		
		// load sequences into an HashMap object
		// TODO memory optimization: load only needed sequences
		SeqlibFile seqlibfile = new SeqlibFile(seqlib);
		HashMap<String, Sequence> library = seqlibfile.getLibrary();
		
		LocalSequenceGotoh gotoh =
				new LocalSequenceGotoh(
					-10.0, -2.0,
					bioinfo.alignment.matrices.QuasarMatrix.DAYHOFF_MATRIX
				);
// Print out Header		
		System.out.println("target\ttemplate\t\"HubeRDP RMSD\"\t\"HubeRDP TMScore\"\t\"HubeRDP GDT\"\t\"Gotoh RMSD\"\t\"Gotoh TMScore\"\t\"Gotoh GDT\"\t\"HubeRDP tree depth\"");
		
		// initialise loop stuff
		String line = "";
		int rdpdepth = 0;
		
		// initialize HubeRDP stuff
		// TODO
		HubeRDP rdp = new HubeRDP();
		rdp.addOracle(new RDPOracle());
		SipplContactPotential sippl = new SipplContactPotential();
		sippl.readFromVPOTFile(vpotfile);
		rdp.setScoring(new RDPScoring(RDPScoring.GAMMA, RDPScoring.DELTA,
				RDPScoring.EPSILON, RDPScoring.ZETA,
				QuasarMatrix.DAYHOFF_MATRIX, new HydrophobicityMatrix(
						hydromatrixFile), new CCPMatrix(ccpfilename), sippl,
				null, RDPScoring.VOROPATH, RDPScoring.GRID_EXTEND,
				RDPScoring.GRID_DENSITY, RDPScoring.GRID_CLASH,
				RDPScoring.MIN_CONTACT));
		
		// initialize TM stuff
		TMMain tmmain = new TMMain();
		
		for (String[] job : joblist) {
			try {
				// load new job
				String templatePDBID = job[0].substring(0, 4);
				String targetPDBID = job[1].substring(0, 4);
				Sequence template = library.get(job[0]);
				Sequence target = library.get(job[1]);
				SequenceAlignment ali =
						new SequenceAlignment(
								target, target,
								target.getSequence(), target.getSequence(),
								0.0
						);
				PDBEntry templateStructure =
					new PDBFileReader(pdbpath).readPDBFromFile(
							PDBFile.getFile(pdbpath, templatePDBID),
							job[0].charAt(4)
					);
				PDBEntry targetStructure =
						new PDBFileReader(pdbpath).readPDBFromFile(
								PDBFile.getFile(pdbpath, targetPDBID),
								job[1].charAt(4)
						);
				
//output				target			template
				line = job[1] + "\t" + job[0]+"\t";
				
			// align sequences with HubeRDP
				// initialize HubeRDP
				System.out.print("Constructing RDP Tree structure...");
				int[][] rows = new int[2][templateStructure.length() + target.length()];
				for(int i = 0; i < template.length(); i++) {
					rows[0][i] = i;
					rows[1][i] = -1;
				}
				for(int i = 0; i < target.length(); i++) {
					rows[0][templateStructure.length()+i] = -1;
					rows[1][+templateStructure.length()+i] = i;
				}
				
				RDPProblem root = new RDPProblem(new Threading(templateStructure, target, rows, 0), 0, rows[0].length-1);
				RDPSolutionTree t = new RDPSolutionTree(root);
				RDPPriorityQueue pq = new RDPPriorityQueue(t.getRoot());
				// run HubeRDP
				rdp.rdp(t, pq);
				// save for output
				rdpdepth = t.getDepth()/2;
				
				// get HubeRDP's alignment
				SequenceAlignment rdpAlignment =
						t.getRoot().getTA().get(0).getThreading().asSequenceAlignment();
			// use CoordMapper on HubeRDP Alignment
				PDBEntry rdpStructure =
						SimpleCoordMapper.map(templateStructure, rdpAlignment);
				
			// calculate TM stuff
				Transformation rdptmtr = tmmain.calculateTransformation(
						ali, rdpStructure, targetStructure);
				
// output
				line += df.format(rdptmtr.getRmsd()) + "\t" +
						df.format(rdptmtr.getTmscore()) + "\t" +
						df.format(rdptmtr.getGdt()) + "\t";
				
			// align Sequences with Gotoh
				SequenceAlignment gotohAlignment =
						gotoh.align(template, target);
				
			// use CoordMapper on Gotoh Alignment
				PDBEntry gotohStructure =
						SimpleCoordMapper.map(templateStructure, gotohAlignment);

				// calculate TM stuff
				Transformation gotohtmtr = tmmain.calculateTransformation(
						ali, gotohStructure, targetStructure);
				
// output
				line += df.format(gotohtmtr.getRmsd()) + "\t" +
						df.format(gotohtmtr.getTmscore()) + "\t" +
						df.format(gotohtmtr.getGdt()) + "\t";
				
				line += rdpdepth;
				System.out.println(line);
				
			} catch (Exception e) {
				System.err.println("Error while working on job "+ job[0]+" " + job[1]);
				e.printStackTrace();
			}
			
		}
		
	}

}

/******************************************************************************
 * "A question that sometimes drives me hazy:                                 *
 *  Am I or are the others crazy?"                                            *
 *     - Albert Einstein (1879 - 1955)                                        *
 ******************************************************************************/
