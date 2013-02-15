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

	// set test data
	private static final Sequence template = new Sequence("1j2xA00","GPLDVQVTEDAVRRYLTRKPMTTKDLLKKFQTKKTGLSSEQTVNVLAQILKRLNPERKMINDKMHFSLK");
	private static final Sequence target   = new Sequence("1dp7P00","TVQWLLDNYETAEGVSLPRSTLYNHYLLHSQEQKLEPVNAASFGKLIRSVFMGLRTRRLGTRGNSKYHYYGLRIK");
//	private static final Sequence template = new Sequence("test001","GGGGCA");
//	private static final Sequence target   = new Sequence("test002","TTTGGGGA");
	
	/**
	 * main function
	 * @param args no args needed
	 */
	public static void main(String[] args) {
		
		// set test data
		PDBEntry templateStructure = null;
		PDBEntry targetStructure = null;
		SequenceAlignment ali = null;
		
		// read test pdb file
		System.out.print("Reading Template structure file...");
		PDBFileReader fr = new PDBFileReader();
		templateStructure = fr.readPDBFromFile
				("/home/h/huberste/gobi/webserver/pdb/1J2XA00.pdb");
		// nullify fr for the Garbage Collector
		fr = null;
		System.out.println(" done!");
		
		// construct rdp tree
		System.out.print("Constructing RDP Tree structure...");
		RDPProblem root = new RDPProblem
				 (template, templateStructure,
					target, targetStructure,
				 ali, 0, template.length() - 1, 0, target.length() - 1);
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
//		rdp.addOracle(new ManualOracle());
		
		// set scoring
		rdp.setScoring(new SimpleScoring());
//		rdp.setScoring(new ManualScoring());
		
		// execute rdp algorithm
		System.out.println("HubeRDP will now be executed!");
		rdp.rdp(t, pq);
		System.out.println("HubeRDP was successfully executed!");
		// Solution is now in t.getRoot();
		
		System.out.println(t.getRoot().getTA().toString());

	}

}

/******************************************************************************
 * "A question that sometimes drives me hazy:                                 *
 *  Am I or are the others crazy?"                                            *
 *     - Albert Einstein (1879 - 1955)                                        *
 ******************************************************************************/
