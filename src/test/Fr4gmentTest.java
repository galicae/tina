package test;

import java.util.ArrayList;
import java.util.LinkedList;
import java.util.List;

import bioinfo.alignment.matrices.QuasarMatrix;
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
		PDBEntry pdb1 = reader.readPDBFromFile("proteins/1BDS.pdb");
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
		String query = "";
		for(int i = 0; i < pdb1.length(); i++) {
			query += pdb1.getAminoAcid(i).getName().getOneLetterCode();
		}
		String curSub = query.substring(0, 5);
		ProteinFragment curFrag = findFragment(curSub, clusters);
		LinkedList<ProteinFragment> result = new LinkedList<ProteinFragment>();
		result.add(curFrag);
		for(int i = 2; i < query.length(); i += 3) {
			curFrag = findFragment(curSub, clusters);
			result.add(curFrag);
		}
	}
	
	public static ProteinFragment findFragment(String query, LinkedList<FragmentCluster> clusters) {
		ProteinFragment curFrag = new ProteinFragment("dada", new double[1][1], 0 ,1);
		double tempScore = Double.MIN_VALUE;
		double temp = 0;
		double[][] matrix = new double[1][1];
		matrix = QuasarMatrix.DAYHOFF_MATRIX;
		
		for(FragmentCluster c: clusters) {
			for(int i = 0; i < c.getCentroid().fragLength; i++) {
				int x = 65 - c.getCentroid().getSequence().charAt(i);
				int y = 65 - query.charAt(i);
				temp += matrix[x][y] * c.getPssm()[x][y];
			}
			if(temp > tempScore) {
				curFrag = c.getCentroid().clone();
				curFrag.setSequence(query);
				tempScore = temp;
			}
		}
		return curFrag;
	}

}
