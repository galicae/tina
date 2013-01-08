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
		// first read everything there is to read
		PDBFileReader reader = new PDBFileReader();
		// and save it in a files list
		List<PDBEntry> files = new LinkedList<PDBEntry>();
		ArrayList<ProteinFragment> pList = new ArrayList<ProteinFragment>();
		PDBEntry pdb1 = reader.readPDBFromFile("1x2tA00.pdb");
		
		files.add(pdb1);
		
		// next make fragments out of all PDBs and save them in pList
		for(PDBEntry e: files) {
			Fragmenter.crunchBackboneSeq(e, pList, 5);
		}
		int initSum = pList.size();
		
		// now the clustering algorithm with its many steps
		DBScan clustah = new DBScan();
		// define cluster list
		LinkedList<FragmentCluster> clusters = new LinkedList<FragmentCluster>();
		clustah.oppaDBStyle(4, 1.0, pList, clusters);
		
		int sumOfFrags = 0;
		for(FragmentCluster c: clusters) {
			sumOfFrags += c.getSize();
		}
		
		// report clustering efficiency
		System.out.format("%d out of %d fragments in %d clusters.\n" , sumOfFrags, initSum, clusters.size());
		
		// write out clusters
		for (FragmentCluster c : clusters) {
			try {
				BufferedWriter br = new BufferedWriter(new FileWriter("/home/p/papadopoulos/Desktop/results/KM/" + c.getCentroid().getID()));
				br.write(c.toString());
				br.close();
			} catch (Exception e) {
				e.printStackTrace();
			}
		}
	}
}
