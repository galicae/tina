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

/**
 * @author huberste
 * @lastchange 2013-01-29
 * 
 */
public class huberdptest {

	private static final Sequence target = new Sequence("1dp7P00","TVQWLLDNYETAEGVSLPRSTLYNHYLLHSQEQKLEPVNAASFGKLIRSVFMGLRTRRLGTRGNSKYHYYGLRIK");
	private static final Sequence template = new Sequence("1j2xA00","GPLDVQVTEDAVRRYLTRKPMTTKDLLKKFQTKKTGLSSEQTVNVLAQILKRLNPERKMINDKMHFSLK");
	
	/**
	 * @param args
	 */
	public static void main(String[] args) {
		
		PDBEntry targetStructure = null;
		PDBEntry templateStructure = null;
		SequenceAlignment ali = null;
		
		PDBFileReader fr = new PDBFileReader();
		templateStructure = fr.readPDBFromFile
				("/home/h/huberste/gobi/webserver/pdb/1J2XA00.pdb");
		RDPProblem root = new RDPProblem
				(target, targetStructure,
				 template, templateStructure,
				 ali, 0, target.length(), 0, template.length());
		
		RDPSolutionTree t = new RDPSolutionTree(root);
		RDPPriorityQueue pq = new RDPPriorityQueue(t.getRoot());
		
//		Oracle oracle = new Oracle();
		
		HubeRDP rdp = new HubeRDP();
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