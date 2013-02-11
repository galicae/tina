package bioinfo.proteins.fragm3nt;

import java.util.LinkedList;

import bioinfo.superpos.*;

import bioinfo.proteins.PDBEntry;

public class CheatAssembler extends Assembler {
	PDBEntry control;

	public CheatAssembler(int i, PDBEntry control) {
		super(i);
		this.control = control;
	}

	public ProteinFragment proveConcept(String query,
			LinkedList<FragmentCluster> clusters, int extent) {

		ProteinFragment resultFragment = assembleProtein(clusters, extent,
				query);
		return resultFragment;
	}

	private ProteinFragment assembleProtein(
			LinkedList<FragmentCluster> clusters, int extent, String query) {

		ProteinFragment tempResult = new ProteinFragment("bla",
				new double[0][0], extent);
		
		ProteinFragment controlFrag = new ProteinFragment("control",
				PDBReduce.reduceSinglePDB(control), extent);
		controlFrag.setSequence(query);
		
		ProteinFragment resultFragment = new ProteinFragment("bla",
				new double[0][0], extent);
		resultFragment = controlFrag.getPart(0, fragLength);
		int size = clusters.size();

		int index = fragLength;
		// 8 5
		int indexLimit = query.length() - fragLength + extent;
		while (index < indexLimit) {
			ProteinFragment temp = clusters.get(0).getCentroid().clone();

			double clusterScore = Double.MAX_VALUE;
			for (int i = 0; i < size; i++) {
				temp = resultFragment.clone();
				clusters.get(i).getCentroid().setSequence(query.substring(index - extent, index - extent + fragLength));
				alignFragments(temp, clusters.get(i).getCentroid(), extent);
				double pairScore = compare(controlFrag.getPart(0, temp.getAllResidues().length), temp);
				if (clusterScore > pairScore) {
					clusterScore = pairScore;
					tempResult = temp.clone();
				}
			}
			resultFragment = tempResult.clone();
			index = resultFragment.getAllResidues().length;
//			System.out.println("index at " + index);
		}
		
		int newExtent = fragLength - (query.length() - resultFragment.getAllResidues().length);
		double clusterScore = Double.MAX_VALUE;
		for (int i = 0; i < size; i++) {
			ProteinFragment temp = clusters.get(0).getCentroid().clone();
			temp = resultFragment.clone();
			clusters.get(i).getCentroid().setSequence(query.substring(index - newExtent, index - newExtent + fragLength));
			alignFragments(temp, clusters.get(i).getCentroid(), newExtent);
			double pairScore = compare(controlFrag.getPart(0, temp.getAllResidues().length), temp);
			if (clusterScore > pairScore) {
				clusterScore = pairScore;
				tempResult = temp.clone();
			}
		}
		resultFragment = tempResult.clone();
		index = resultFragment.getAllResidues().length;
//		System.out.println("index at " + index);
		
		resultFragment.setSequence(query);
		return resultFragment;
	}

	private double compare(ProteinFragment nativStr, ProteinFragment predStr) {
		double[][][] kabschFood = new double[2][predStr.getAllResidues().length][3];
		kabschFood[0] = nativStr.getAllResidues();
		kabschFood[1] = predStr.getAllResidues();
		Transformation t = Kabsch.calculateTransformation(kabschFood);
		return t.getRmsd();
	}

}
