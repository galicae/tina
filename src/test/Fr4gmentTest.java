package test;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.util.LinkedList;
import java.util.List;

import bioinfo.alignment.matrices.QuasarMatrix;
import bioinfo.proteins.PDBEntry;
import bioinfo.proteins.PDBFileReader;
import bioinfo.proteins.fragm3nt.FragmentCluster;
import bioinfo.proteins.fragm3nt.Fragmenter;
import bioinfo.proteins.fragm3nt.KMeansAllvsAll;
import bioinfo.proteins.fragm3nt.ProteinFragment;
import bioinfo.superpos.Kabsch;
import bioinfo.superpos.Transformation;

public class Fr4gmentTest {
	static int fragLength = 5;
	public static void main(String[] args) {
		PDBFileReader reader = new PDBFileReader("./proteins/");
		List<PDBEntry> files = new LinkedList<PDBEntry>();
		LinkedList<ProteinFragment> pList = new LinkedList<ProteinFragment>();
		files = reader.readPdbFolder();
		// files.add(pdb1);
		PDBEntry pdb1 = files.get(0);
		for (PDBEntry e : files) {
			Fragmenter.crunchBackboneSeq(e, pList, 5);
		}
		KMeansAllvsAll clustah = new KMeansAllvsAll(pList);
		LinkedList<FragmentCluster> clusters = new LinkedList<FragmentCluster>();
		clustah.initializeClusters();
		clustah.update(100);

		clusters = clustah.getClusters();
		
		for(FragmentCluster fc: clusters) {
			for(ProteinFragment p: fc.getFragments()) {
				p.correctCoordinates();
			}
		}

		LinkedList<ProteinFragment> curFrags = new LinkedList<ProteinFragment>();
		double[][] pssm = new double[fragLength][26];
		char c = 'a';
		for (FragmentCluster fr : clusters) {
			pssm = new double[fragLength][26];
			curFrags = fr.getFragments();
			for (ProteinFragment f : curFrags) {
				for (int i = 0; i < f.getSequence().length(); i++) {
					c = f.getSequence().charAt(i);
					pssm[i][c - 65] += (1.0 / fr.getSize());
				}
			}
			fr.setPssm(pssm);
		}
		
		// assembly?
		String query = "";
		for (int i = 0; i < pdb1.length(); i++) {
			query += pdb1.getAminoAcid(i).getName().getOneLetterCode();
		}
		
		String curSub = query.substring(0, 5);
		ProteinFragment curFrag = findFragment(curSub, clusters);
		LinkedList<ProteinFragment> result = new LinkedList<ProteinFragment>();
		result.add(curFrag);
		for (int i = 3; i < query.length() - 5; i += 3) {
			curSub = query.substring(i, i + 5);
			curFrag = findFragment(curSub, clusters);
			result.add(curFrag);
		}
		
		pssm = null;
		curSub = null;
		curFrag = null;
		
		StringBuilder sb = new StringBuilder();
		sb.append(result.get(0).toString());
//		sb.append(result.get(1).toString());
		String temp = "";
		for(int i = 1; i < result.size(); i++) {
			temp = alignFragments(result.get(i - 1), result.get(i), 2, 5 + 3*i);
			sb.append(temp.toString());
		}
		
		try {
			BufferedWriter w = new BufferedWriter(new FileWriter("9pap_test.pdb"));
			w.write(sb.toString());
			w.close();
		}
		catch(Exception e) {
			e.printStackTrace();
		}
	}

	public static ProteinFragment findFragment(String query,
			LinkedList<FragmentCluster> clusters) {
		ProteinFragment curFrag = new ProteinFragment("dada", new double[1][1],
				0, 1);
		double tempScore = - Double.MAX_VALUE;
		double temp = 0;
		double[][] matrix = new double[1][1];
		matrix = QuasarMatrix.DAYHOFF_MATRIX;
		
		for (FragmentCluster c : clusters) {
			temp = 0;
			for (int i = 0; i < fragLength / 4; i++) {
				int y = query.charAt(i) - 65;
				temp += matrix[i][y] * c.getPssm()[i][y];
			}
			if (temp > tempScore) {
				curFrag = c.getCentroid().clone();
				curFrag.setSequence(query);
				tempScore = temp;
			}
		}
		return curFrag;
	}
	
	public static String alignFragments(ProteinFragment stable, ProteinFragment move, int extent, int position) {
		double[][][] kabschFood = new double[2][extent * 4][3];
		int shove = (fragLength - extent) * 4;
		for (int i = 0; i < extent * 4; i++) {
			kabschFood[0][i] = stable.getResidue(shove + i);
			kabschFood[1][i] = move.getResidue(i);
		}
		Transformation t = Kabsch.calculateTransformation(kabschFood);
		move.setCoordinates(t.transform(move.getAllResidues()));
		matchCoordinates(stable, move, extent, position);
		return move.toString(extent, position);
	}

	private static void matchCoordinates(ProteinFragment stable,
			ProteinFragment move, int extent, int position) {
		int last = stable.getAllResidues().length - 1;
		double[] correct = new double[3];
		double[] lastStable = stable.getAllResidues()[last];
		double[] firstMove = move.getAllResidues()[extent * 4 - 1];
		for(int i = 0; i < 3; i++) {
			correct[i] = lastStable[i] - firstMove[i];
		}
//		firstMove = move.getAllResidues()[extent * 4 - 1];
//		System.err.println(firstMove[0] + "\t" + firstMove[1] + "\t" + firstMove[2] + "\t");
		move.translateCoordinates(correct);
//		System.err.println(firstMove[0] + "\t" + firstMove[1] + "\t" + firstMove[2] + "\t");
	}

}
