package test;

import java.util.LinkedList;

import bioinfo.pdb.PDBFile;
import bioinfo.proteins.PDBEntry;
import bioinfo.proteins.PDBFileReader;
import bioinfo.proteins.framg3nt.Fragmenter;
import bioinfo.proteins.framg3nt.KMeansAllvsAll;
import bioinfo.proteins.framg3nt.ProteinFragment;

public class Fragm3ntTest {
	public static void main(String[] args) {
		PDBFileReader reader = new PDBFileReader();
		PDBEntry pdb1 = reader.readPDBFromFile("1TIMA00.pdb");

		LinkedList<ProteinFragment> pList = new LinkedList<ProteinFragment>();
		Fragmenter.crunch(pdb1, pList, 5);
		
		KMeansAllvsAll clustah = new KMeansAllvsAll(pList);
		clustah.initializeClusters();
		System.out.println("initialized clusters");
		clustah.toTextFiles("init");
		clustah.update();
		System.out.println("updated");
		clustah.toTextFiles("upd");
	}
}
