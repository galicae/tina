package test;

import java.util.LinkedList;

import bioinfo.proteins.PDBEntry;
import bioinfo.proteins.PDBFileReader;
import bioinfo.proteins.framg3nt.FragmentCluster;
import bioinfo.proteins.framg3nt.Fragmenter;
import bioinfo.proteins.framg3nt.KMeansAllvsAll;
import bioinfo.proteins.framg3nt.ProteinFragment;

public class Fragm3ntTest {
	public static void main(String[] args) {
		PDBFileReader reader = new PDBFileReader();
		PDBEntry pdb1 = reader.readPDBFromFile("1TIMA00.pdb");

		LinkedList<ProteinFragment> pList = new LinkedList<ProteinFragment>();
		Fragmenter.crunchBackboneN(pdb1, pList, 7);
		int initSum = pList.size();
		
		KMeansAllvsAll clustah = new KMeansAllvsAll(pList);
		clustah.initializeClusters();
		System.out.println("initialized clusters");
//		clustah.toTextFiles("init");
		
		int sumOfFrags = 0;
		for(FragmentCluster c: clustah.getClusters()) {
			sumOfFrags += c.getSize();
		}
		System.out.format("%d out of %d fragments in %d clusters.\n" , sumOfFrags, initSum, clustah.getClusters().size());
		clustah.update(20);
		sumOfFrags = 0;
		for(FragmentCluster c: clustah.getClusters()) {
			sumOfFrags += c.getSize();
		}
		System.out.println("updated");
		clustah.toTextFiles("upd");
		System.err.format("%d out of %d fragments in %d clusters.\n" , sumOfFrags, initSum, clustah.getClusters().size());
	}
}
