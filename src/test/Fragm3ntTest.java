package test;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.util.ArrayList;
import java.util.LinkedList;
import java.util.List;

import bioinfo.proteins.PDBEntry;
import bioinfo.proteins.PDBFileReader;
import bioinfo.proteins.fragm3nt.DBScan;
import bioinfo.proteins.fragm3nt.FragmentCluster;
import bioinfo.proteins.fragm3nt.Fragmenter;
import bioinfo.proteins.fragm3nt.ProteinFragment;
import bioinfo.superpos.Kabsch;
import bioinfo.superpos.Transformation;

public class Fragm3ntTest {
	public static void main(String[] args) {
//		try {
//			BufferedReader br = new BufferedReader(new FileReader("pdb_catalogue.txt"));
//			String line = "";
//			while((line = br.readLine()) != null) {
//				PDBFile.downloadPDB(line, "./proteins/");
//				System.out.println("got " + line);
//			}
//		}
//		catch(Exception e) {
//			e.printStackTrace();
//		}
		
		PDBFileReader reader = new PDBFileReader();
		
		List<PDBEntry> files = new LinkedList<PDBEntry>();
		ArrayList<ProteinFragment> pList = new ArrayList<ProteinFragment>();
		PDBEntry pdb1 = reader.readPDBFromFile("1x2tA00.pdb");
		
		files.add(pdb1);
		for(PDBEntry e: files) {
			Fragmenter.crunchBackboneSeq(e, pList, 5);
		}
		int initSum = pList.size();
		
		DBScan clustah = new DBScan();
		LinkedList<FragmentCluster> clusters = new LinkedList<FragmentCluster>();
		clustah.oppaDBStyle(4, 1.0, pList, clusters);
//		clustah.toTextFiles("init");
		
		int sumOfFrags = 0;
		for(FragmentCluster c: clusters) {
			sumOfFrags += c.getSize();
		}
		System.out.format("%d out of %d fragments in %d clusters.\n" , sumOfFrags, initSum, clusters.size());
//		clustah.update(20);
//		sumOfFrags = 0;
//		for(FragmentCluster c: clustah.getClusters()) {
//			sumOfFrags += c.getSize();
//		}
		
		for (FragmentCluster c : clusters) {
			try {
				BufferedWriter br = new BufferedWriter(new FileWriter("DBS_" + c.getCentroid().getID()));
				br.write(c.toString());
				br.close();
			} catch (Exception e) {
				e.printStackTrace();
			}
		}
	}
}
