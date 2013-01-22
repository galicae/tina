package bioinfo.proteins.fragm3nt;

import java.util.LinkedList;
import java.util.List;

import bioinfo.alignment.matrices.QuasarMatrix;
import bioinfo.proteins.PDBEntry;
import bioinfo.superpos.Kabsch;
import bioinfo.superpos.Transformation;

public class Assembler {
	private final int fragLength;
	
	public Assembler(int fragLength) {
		this.fragLength = fragLength;
	}

	private ProteinFragment findFragment(String query,
			LinkedList<FragmentCluster> clusters) {
		ProteinFragment curFrag = new ProteinFragment("dada", new double[1][1], 1);
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
	
	private void alignFragments(ProteinFragment stable, ProteinFragment move, int extent) {
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

	private void matchCoordinates(ProteinFragment stable, ProteinFragment move, int extent) {
		int last = stable.getAllResidues().length - 1;
		double[] correct = new double[3];
		double[] lastStable = stable.getAllResidues()[last];
		double[] firstMove = move.getAllResidues()[extent - 1];
		for(int i = 0; i < 3; i++) {
			correct[i] = lastStable[i] - firstMove[i];
		}
		move.translateCoordinates(correct);
	}
	
	public String readSequence(PDBEntry pdb1) {
		String query = "";
		for (int i = 0; i < pdb1.length(); i++) {
			query += pdb1.getAminoAcid(i).getName().getOneLetterCode();
		}
		return query;
	}
	
	private LinkedList<ProteinFragment> collectFragments(String query, LinkedList<FragmentCluster> clusters) {
		String curSub = query.substring(0, 5);
		ProteinFragment curFrag = findFragment(curSub, clusters);
		LinkedList<ProteinFragment> result = new LinkedList<ProteinFragment>();
		result.add(curFrag);
		for (int i = 3; i < query.length() - 5; i += 3) {
			curSub = query.substring(i, i + 5);
			curFrag = findFragment(curSub, clusters);
			result.add(curFrag);
		}
		return result;
	}
	
	private ProteinFragment assembleProtein(List<ProteinFragment> rightFragments) {
		ProteinFragment resultFragment = new ProteinFragment("res", new double[1][1], fragLength);
		resultFragment = rightFragments.get(0).clone();
		for(int i = 1; i < rightFragments.size(); i++) {
			alignFragments(resultFragment, rightFragments.get(i), 2);
		}
		return resultFragment;
	}
	
	public ProteinFragment predictStructure(String query, LinkedList<FragmentCluster> clusters) {
		LinkedList<ProteinFragment> result = new LinkedList<ProteinFragment>();
		result = collectFragments(query, clusters);
		
		ProteinFragment resultFragment = assembleProtein(result);
		return resultFragment;
	}
}
