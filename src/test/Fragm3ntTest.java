package test;

import java.io.BufferedReader;
import java.io.FileReader;

import bioinfo.pdb.PDBFile;

public class Fragm3ntTest {
	public static void main(String[] args) {
		try {
			BufferedReader br = new BufferedReader(new FileReader("pdb_catalogue.txt"));
			String line = "";
			while((line = br.readLine()) != null) {
				PDBFile.downloadPDB(line, "./proteins/");
				System.out.println("got " + line);
			}
		}
		catch(Exception e) {
			e.printStackTrace();
		}
		
//		PDBFileReader reader = new PDBFileReader("./proteins/");
//		
//		List<PDBEntry> files = new LinkedList<PDBEntry>();
//		LinkedList<ProteinFragment> pList = new LinkedList<ProteinFragment>();
//		files = reader.readPdbFolder();
//		for(PDBEntry e: files) {
//			Fragmenter.crunchBackboneN(e, pList, 7);
//		}
//		int initSum = pList.size();
//		
//		KMeansAllvsAll clustah = new KMeansAllvsAll(pList);
//		clustah.initializeClusters();
//		System.out.println("initialized clusters");
////		clustah.toTextFiles("init");
//		
//		int sumOfFrags = 0;
//		for(FragmentCluster c: clustah.getClusters()) {
//			sumOfFrags += c.getSize();
//		}
//		System.out.format("%d out of %d fragments in %d clusters.\n" , sumOfFrags, initSum, clustah.getClusters().size());
//		clustah.update(2000);
//		sumOfFrags = 0;
//		for(FragmentCluster c: clustah.getClusters()) {
//			sumOfFrags += c.getSize();
//		}
//		System.out.println("updated");
//		clustah.toTextFiles("upd");
//		System.err.format("%d out of %d fragments in %d clusters.\n" , sumOfFrags, initSum, clustah.getClusters().size());
	}
}
