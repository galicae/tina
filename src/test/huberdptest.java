/******************************************************************************
 * hubeRDPTest is mainly a test routine for calling hubeRDP.                  *
 *                                                                            *
 * This file is best read at line width 80 and tab width 4.                   *
 *                                                                   huberste *
 ******************************************************************************/
package test;

import bioinfo.Sequence;
import bioinfo.alignment.SequenceAlignment;
import bioinfo.proteins.PDBEntry;
import bioinfo.proteins.PDBFileReader;
import huberdp.*;
import huberdp.oracles.TinyOracle;

/**
 * @author huberste
 * @lastchange 2013-01-29
 * 
 */
public class huberdptest {

	// set test data
	private static final Sequence target = new Sequence("1dp7P00","TVQWLLDNYETAEGVSLPRSTLYNHYLLHSQEQKLEPVNAASFGKLIRSVFMGLRTRRLGTRGNSKYHYYGLRIK");
	private static final Sequence template = new Sequence("1j2xA00","GPLDVQVTEDAVRRYLTRKPMTTKDLLKKFQTKKTGLSSEQTVNVLAQILKRLNPERKMINDKMHFSLK");
	
	/**
	 * main function
	 * @param args
	 */
	public static void main(String[] args) {
		
		// set test data
		PDBEntry targetStructure = null;
		PDBEntry templateStructure = null;
		SequenceAlignment ali = null;
		
		// read test pdb file
		System.out.print("Reading Template structure file...");
		PDBFileReader fr = new PDBFileReader();
		templateStructure = fr.readPDBFromFile
				("/home/h/huberste/gobi/webserver/pdb/1J2XA00.pdb");
		System.out.println(" done!");
		
		// construct rdp tree
		System.out.print("Constructing RDP Tree structure...");
		RDPProblem root = new RDPProblem
				(target, targetStructure,
				 template, templateStructure,
				 ali, 0, target.length(), 0, template.length());
		RDPSolutionTree t = new RDPSolutionTree(root);
		System.out.println(" done!");

		// construct priority queue
		System.out.print("constructing Priority Queue...");
		RDPPriorityQueue pq = new RDPPriorityQueue(t.getRoot());
		System.out.println(" done!");
		
		// construct RDP
		HubeRDP rdp = new HubeRDP();
		
		// add oracles
		rdp.addOracle(new TinyOracle());
		
		// execute rdp algorithm
		System.out.println("HubeRDP will now be executed!");
		rdp.rdp(t, pq);
		System.out.println("HubeRDP was successfully executed!");
		// Solution is now in t.getRoot();
		
		// TODO output of t.getRoot();

	}

}

/******************************************************************************
 * "A question that sometimes drives me hazy:                                 *
 *  Am I or are the others crazy?"                                            *
 *     - Albert Einstein (1979 - 1966)                                        *
 ******************************************************************************/