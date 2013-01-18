package test;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.util.LinkedList;
import java.util.List;

import bioinfo.alignment.matrices.QuasarMatrix;
import bioinfo.proteins.Atom;
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
		System.out.println("crunched");
		KMeansAllvsAll clustah = new KMeansAllvsAll(pList);
		LinkedList<FragmentCluster> clusters = new LinkedList<FragmentCluster>();
		clustah.initializeClusters();
		clustah.update(200);

		clusters = clustah.getClusters();
		
//		for(FragmentCluster fc: clusters) {
//			for(ProteinFragment p: fc.getFragments()) {
//				p.correctCoordinates();
//			}
//		}

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
		
		ProteinFragment resultFragment = new ProteinFragment("res", new double[1][1], 0, fragLength);
		resultFragment = result.get(0).clone();
		
		for(int i = 1; i < 2; i++) {
			alignFragments(resultFragment, result.get(i), 2);
//			resultFragment.append(result.get(i).getAllResidues(), result.get(i).getSequence());
		}
		
		try {
			BufferedWriter w = new BufferedWriter(new FileWriter("C:/users/nikos/Desktop/test.pdb"));
			w.write(resultFragment.toString());
			w.close();
		}
		catch(Exception e) {
			e.printStackTrace();
		}
	}

	public static ProteinFragment findFragment(String query,
			LinkedList<FragmentCluster> clusters) {
		ProteinFragment curFrag = new ProteinFragment("dada", new double[1][1], 0, 1);
		double tempScore = - Double.MAX_VALUE;
		double temp = 0;
		double[][] matrix = new double[1][1];
		matrix = QuasarMatrix.DAYHOFF_MATRIX;
		
		for (FragmentCluster c : clusters) {
			temp = 0;
			for (int i = 0; i < fragLength; i++) {
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
	
	public static void alignFragments(ProteinFragment stable, ProteinFragment move, int extent) {
		double[][][] kabschFood = new double[2][extent][3];
		int shove = stable.getAllResidues().length - extent;
		for (int i = 0; i < extent; i++) {
			kabschFood[0][i] = stable.getResidue(shove + i);
			kabschFood[1][i] = move.getResidue(i);
		}
		Transformation t = Kabsch.calculateTransformation(kabschFood);
		move.setCoordinates(t.transform(move.getAllResidues()));
		matchCoordinates(stable, move, extent);
		
		double[][] moveCoordinates = new double[(fragLength - extent)][3];
		double[][] allResidues = move.getAllResidues();
		for(int i = 0; i < moveCoordinates.length; i++) {
			moveCoordinates[i][0] = allResidues[extent + i][0];
			moveCoordinates[i][1] = allResidues[extent + i][1];
			moveCoordinates[i][2] = allResidues[extent + i][2];
		}
		String moveSeq = move.getSequence().substring(extent);
		
		stable.append(moveCoordinates, moveSeq);
	}

	private static void matchCoordinates(ProteinFragment stable, ProteinFragment move, int extent) {
		int last = stable.getAllResidues().length - 1;
		double[] correct = new double[3];
		double[] lastStable = stable.getAllResidues()[last];
		double[] firstMove = move.getAllResidues()[extent - 1];
		for(int i = 0; i < 3; i++) {
			correct[i] = lastStable[i] - firstMove[i];
		}
//		firstMove = move.getAllResidues()[extent - 1];
//		System.err.println(firstMove[0] + "\t" + firstMove[1] + "\t" + firstMove[2] + "\t");
		move.translateCoordinates(correct);
//		System.err.println(firstMove[0] + "\t" + firstMove[1] + "\t" + firstMove[2] + "\t");
	}

}
