package test;

import java.util.ArrayList;
import java.util.LinkedList;
import java.util.List;

import bioinfo.proteins.PDBEntry;
import bioinfo.proteins.PDBFileReader;
import bioinfo.proteins.fragm3nt.DBScan;
import bioinfo.proteins.fragm3nt.FragmentCluster;
import bioinfo.proteins.fragm3nt.Fragmenter;
import bioinfo.proteins.fragm3nt.KMeansAllvsAll;
import bioinfo.proteins.fragm3nt.ProteinFragment;

public class Fr4gmentTest {

	public static void main(String[] args) {
		PDBFileReader reader = new PDBFileReader();
		List<PDBEntry> files = new LinkedList<PDBEntry>();
		LinkedList<ProteinFragment> pList = new LinkedList<ProteinFragment>();
		PDBEntry pdb1 = reader.readPDBFromFile("1x2tA00.pdb");
		files.add(pdb1);
		for (PDBEntry e : files) {
			Fragmenter.crunchBackboneSeq(e, pList, 5);
		}
		KMeansAllvsAll clustah = new KMeansAllvsAll(pList);
		LinkedList<FragmentCluster> clusters = new LinkedList<FragmentCluster>();
		clustah.initializeClusters();
		clustah.update(200);
		
		clusters = clustah.getClusters();
		
		LinkedList<ProteinFragment> curFrags = new LinkedList<ProteinFragment>();
		int fragLength = clusters.getFirst().getFragments().getFirst().fragLength;
		double[][] pssm = new double[fragLength][26];
		char c = 'a';
		for(FragmentCluster fr: clusters) {
			pssm = new double[fragLength][26];
			curFrags = fr.getFragments();
			for(ProteinFragment f: curFrags) {
				for(int i = 0; i < f.getSequence().length(); i++) {
					c = f.getSequence().charAt(i);
					pssm[i][c-65] += (1.0/fr.getSize());
				}
			}
			fr.setPssm(pssm);
		}
		
		// assembly?
		String query = "DCPSGWSSYEGHCYKPFKLYKTWDDAERFCTEQAKGGHLVSIESAGEADFVAQLVTENIQNTKSYVWIGLRVQGKEKQCSSEWSDGSSVSYENWIEAESKTCLGLEKETGFRKWVNIYCGQQNPFVCEA";
		String curSub = query.substring(0, 5);
		
		for(int i = 2; i < query.length(); i += 3) {
			
		}
	}
	
	public ProteinFragment findFragment(String query, FragmentCluster cluster) {
		
		return null;
	}

}
