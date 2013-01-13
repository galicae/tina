package test;

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
		clustah.update(200);

		clusters = clustah.getClusters();

		LinkedList<ProteinFragment> curFrags = new LinkedList<ProteinFragment>();
		int fragLength = clusters.getFirst().getFragments().getFirst().fragLength;
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

		StringBuilder sb = new StringBuilder();
		sb.append(result.getFirst().toString());
		String temp = "";
		for (int i = 0; i < result.size() - 1; i++) {
			temp = alignFragments(result.get(i), result.get(i + 1), i, fragLength);
			sb.append(temp);
		}
		System.out.println(sb.toString());
	}

	private static String alignFragments(ProteinFragment stable,
			ProteinFragment move, int position, int fragLength) {
		double[][][] kabschFood = new double[2][2 * 4][3];
		for(int i = 0; i < kabschFood[0].length; i++) {
			kabschFood[0][i] = stable.getResidue(i);
			kabschFood[1][i] = move.getResidue(fragLength - 2 + i);
		}
		
//		kabschFood[1] = move.getAllResidues();

		Transformation t = Kabsch.calculateTransformation(kabschFood);
		kabschFood[1] = t.transform(kabschFood[1]);
		move.setCoordinates(kabschFood[1]);
		return move.toString(1, position);
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
			for (int i = 0; i < c.getCentroid().fragLength / 4; i++) {
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

}
