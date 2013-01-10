package test;

import java.util.ArrayList;
import java.util.LinkedList;
import java.util.List;

import bioinfo.proteins.PDBEntry;
import bioinfo.proteins.PDBFileReader;
import bioinfo.proteins.fragm3nt.DBScan;
import bioinfo.proteins.fragm3nt.FragmentCluster;
import bioinfo.proteins.fragm3nt.Fragmenter;
import bioinfo.proteins.fragm3nt.ProteinFragment;

public class Fr4gmentTest {

	public static void main(String[] args) {
		PDBFileReader reader = new PDBFileReader();
		List<PDBEntry> files = new LinkedList<PDBEntry>();
		ArrayList<ProteinFragment> pList = new ArrayList<ProteinFragment>();
		PDBEntry pdb1 = reader.readPDBFromFile("1x2tA00.pdb");
		files.add(pdb1);
		for (PDBEntry e : files) {
			Fragmenter.crunchBackboneSeq(e, pList, 5);
		}
		DBScan clustah = new DBScan();
		LinkedList<FragmentCluster> clusters = new LinkedList<FragmentCluster>();
		clustah.oppaDBStyle(4, 1.0, pList, clusters);
		
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

	}

}
