/******************************************************************************
 * HubeRDPValidator.java                                                      *
 *                                                                            *
 * This file is best read at line width 80 and tab width 4.                   *
 *                                                                   huberste *
 ******************************************************************************/
package validation.huberdp;

import huberdp.HubeRDP;

import java.util.HashMap;
import java.util.LinkedList;
//import java.util.Locale;

import bioinfo.Sequence;
import bioinfo.alignment.SequenceAlignment;
import bioinfo.alignment.gotoh.FreeshiftSequenceGotoh;
import bioinfo.pdb.PDBFile;
import bioinfo.proteins.PDBEntry;
import bioinfo.proteins.PDBFileReader;
import bioinfo.proteins.structure.SimpleCoordMapper;

import files.PairFile;
import files.SeqlibFile;
import util.CommandLineParser;

/**
 * HubeRDPValidator validates HubeRDP as an Sequence Alignment Tool against
 * Gotoh.
 * @author huberste
 * @lastchange 2013-02-12
 */
public class HubeRDPValidator {

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
	"java HubeRDPValidator --seqlib domains.seqlib --pairs cathscop.inpairs";
	
	/**
	 * Validates the HubeRDP.
	 * @param args
	 */
	public static void main(String[] args) {
		
//		Locale.setDefault(Locale.US);
		
		// get command line arguments
		CommandLineParser clp = new CommandLineParser(args);
		String seqlib = clp.getStringArg("--seqlib");
		String pairs = clp.getStringArg("--pairs");
		String pdbpath = clp.getStringArg("--pdbpath");
		/*String matrixname = clp.getStringArg("-m");
		double gapopen = clp.getDoubleArg("-go",-10.0);
		double gapextend = clp.getDoubleArg("-ge", -1.0);
		String mode = clp.getStringArg("-mode", "freeshift");
		boolean printali = clp.getBoolArg("-printali", false);
		String printmatrices = clp.getStringArg("-printmatrices", null);
		boolean check = clp.getBoolArg("-check", false);*/
		
		if (pairs == null || seqlib == null) {
			System.out.println(usage);
			System.exit(1);
		}
		
		// load joblist
		PairFile pairfile = new PairFile(pairs);
		LinkedList<String[]> joblist = pairfile.getJoblist();
		
		// load sequences into an HashMap object
		// TODO memory optimization: load only needed sequences
		SeqlibFile seqlibfile = new SeqlibFile(seqlib);
		HashMap<String, Sequence> library = seqlibfile.getLibrary();
		
		FreeshiftSequenceGotoh gotoh =
				new FreeshiftSequenceGotoh(
					10.0, -2.0,
					bioinfo.alignment.matrices.QuasarMatrix.DAYHOFF_MATRIX
				);
		
		for (String[] job : joblist) {
			String pdbid = job[0].substring(0, 4);
			Sequence template = library.get(job[0]);
			Sequence target = library.get(job[1]);
			PDBEntry templateStructure =
				new PDBFileReader(pdbpath).readPDBFromFile(
					PDBFile.getFile(pdbpath, pdbid), job[0].charAt(4)
				);
			System.out.println(">" + job[0] + " " + job[1]);
			
			// align Sequences with HubeRDP
			System.out.println("HubeRDP:");
			SequenceAlignment rdpAlignment =
					HubeRDP.hubeRDPAlign(template, target);
			System.out.println(rdpAlignment.toStringVerbose());
			
			// use CoordMapper on HubeRDP Alignment
			PDBEntry rdpStructure =
					SimpleCoordMapper.map(rdpAlignment, templateStructure);
			
			// align Sequences with Gotoh
			System.out.println("Gotoh:");
			SequenceAlignment gotohAlignment =
					gotoh.align(template, target);
			System.out.println(gotohAlignment.toStringVerbose());
			
			// use CoordMapper on Gotoh Alignment
			PDBEntry gotohStructure =
					SimpleCoordMapper.map(gotohAlignment, templateStructure);
			
			// TODO measure TMScore for each Alignment
			
			// TODO compare TMScores
			
		}
		
	}

}

/******************************************************************************
 * "A question that sometimes drives me hazy:                                 *
 *  Am I or are the others crazy?"                                            *
 *     - Albert Einstein (1979 - 1966)                                        *
 ******************************************************************************/