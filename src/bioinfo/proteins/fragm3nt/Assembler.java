package bioinfo.proteins.fragm3nt;

import java.util.LinkedList;
import java.util.List;

import bioinfo.alignment.matrices.QuasarMatrix;
import bioinfo.proteins.PDBEntry;
import bioinfo.superpos.Kabsch;
import bioinfo.superpos.Transformation;

public class Assembler {
	private final int fragLength;

	/**
	 * the constructor of an assembler object, that with the appropriate input
	 * will be able to predict the backbone of a protein
	 * 
	 * @param fragLength
	 *            the fragment length of the library used
	 */
	public Assembler(int fragLength) {
		this.fragLength = fragLength;
	}

	/**
	 * this function goes through the library with a query of length
	 * fragmentLength and searches the fragment whose sequence profile best
	 * matches this of the query. This fragment is returned as the best
	 * structure prediction for the sequence fragment. More sophisticated
	 * versions of this function should follow.
	 * 
	 * @param query
	 *            the sequence fragment (.length= fragLength) for which we are
	 *            searching a good entrance in the library
	 * @param clusters
	 *            the aforementioned fragment library
	 * @return a clone of the centroid of the best cluster match to query
	 */
	private ProteinFragment findFragment(String query,
			LinkedList<FragmentCluster> clusters) {
		ProteinFragment curFrag = new ProteinFragment("dada", new double[1][1],
				1);
		double tempScore = -Double.MAX_VALUE;
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

	/**
	 * this function merges two fragments by aligning the end of the stable
	 * fragment to the beginning of the "move" fragment and incorporating the
	 * unaligned parts of the "move" fragment to the stable.
	 * 
	 * @param stable
	 *            the fragment we won't move
	 * @param move
	 *            the fragment that is to be moved on the stable. It will be
	 *            optimally rotated so that its beginning will coincide with the
	 *            end of the stable fragment.
	 * @param extent
	 *            how many points are to be aligned
	 */
	private void alignFragments(ProteinFragment stable, ProteinFragment move,
			int extent) {
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
		for (int i = 0; i < moveCoordinates.length; i++) {
			moveCoordinates[i][0] = allResidues[extent + i][0];
			moveCoordinates[i][1] = allResidues[extent + i][1];
			moveCoordinates[i][2] = allResidues[extent + i][2];
		}
		String moveSeq = move.getSequence().substring(extent);

		stable.append(moveCoordinates, moveSeq);
	}

	/**
	 * this function translates the "move" fragment, so that after the rotation
	 * the last aligned residue (the one at {@value extent}) coincides perfectly
	 * with the last residue of the stable fragment; this ensures that the
	 * "move" fragment keeps normal CA distances from the stable
	 * 
	 * @param stable
	 *            the fragment we won't move
	 * @param move
	 *            the fragment that is to be moved on the stable. It will be
	 *            optimally rotated so that its beginning will coincide with the
	 *            end of the stable fragment.
	 * @param extent
	 *            how many points have been aligned
	 */
	private void matchCoordinates(ProteinFragment stable, ProteinFragment move,
			int extent) {
		int last = stable.getAllResidues().length - 1;
		double[] correct = new double[3];
		double[] lastStable = stable.getAllResidues()[last];
		double[] firstMove = move.getAllResidues()[extent - 1];
		for (int i = 0; i < 3; i++) {
			correct[i] = lastStable[i] - firstMove[i];
		}
		move.translateCoordinates(correct);
	}

	/**
	 * helper function for testing purposes. Reads the one-letter amino acid
	 * sequence from a PDBEntry.
	 * 
	 * @param pdb
	 *            the PDBEntry to read
	 * @return the amino acid sequence of the protein
	 */
	public String readSequence(PDBEntry pdb) {
		String query = "";
		for (int i = 0; i < pdb.length(); i++) {
			query += pdb.getAminoAcid(i).getName().getOneLetterCode();
		}
		return query;
	}

	/**
	 * this function goes through the library and collects all fragments that
	 * will be necessary for the assembly of the protein. It does this by
	 * locating the best structure fragments for sequence fragments of length
	 * fragLength, ensuring that the query fragments overlap in {@value extent}
	 * positions, so that the structures can overlap at assembly as well.
	 * 
	 * @param query
	 *            the sequence of the protein whose structure is to be assembled
	 * @param clusters
	 *            the fragment library
	 * @param extent
	 *            how many points have to overlap
	 * @return a list of the correct fragments, that only need to be brought in
	 *         the right relative coordination now
	 */
	private LinkedList<ProteinFragment> collectFragments(String query,
			LinkedList<FragmentCluster> clusters, int extent) {
		String curSub = query.substring(0, fragLength);
		ProteinFragment curFrag = findFragment(curSub, clusters);
		LinkedList<ProteinFragment> result = new LinkedList<ProteinFragment>();
		result.add(curFrag);
		int add = fragLength - extent;
		for (int i = 3; i < query.length() - fragLength; i += add) {
			curSub = query.substring(i, i + fragLength);
			curFrag = findFragment(curSub, clusters);
			result.add(curFrag);
		}
		return result;
	}

	/**
	 * this function assembles a protein, given a list of fragments and the
	 * extent to which they overlap. It does so by aligning the overlapping
	 * regions of two consecutive fragments and adding the non-overlapping
	 * regions to the result.
	 * 
	 * @param rightFragments
	 *            a list of the fragments the algorithm found best for the
	 *            protein sequence
	 * @param extent
	 *            the extent in which points of the fragments overlap
	 * @return the protein structure in the form of a (very large) protein
	 *         fragment.
	 */
	private ProteinFragment assembleProtein(
			List<ProteinFragment> rightFragments, int extent) {
		ProteinFragment resultFragment = new ProteinFragment("res",
				new double[1][1], fragLength);
		resultFragment = rightFragments.get(0).clone();
		for (int i = 1; i < rightFragments.size(); i++) {
			alignFragments(resultFragment, rightFragments.get(i), extent);
		}
		return resultFragment;
	}

	/**
	 * this function summarizes the whole Assembler pipeline. It receives a
	 * query string, a fragment library and the extent to which fragments are
	 * allowed to overlap and returns a protein model.
	 * 
	 * @param query
	 *            the sequence of the predicted protein
	 * @param clusters
	 *            the fragment library in the form generated by the Fragm3nt
	 *            package pipeline
	 * @param extent
	 *            the extent to which fragments are supposed to overlap
	 * @return the protein structure in the form of a (very large) protein
	 *         fragment.
	 */
	public ProteinFragment predictStructure(String query,
			LinkedList<FragmentCluster> clusters, int extent) {
		LinkedList<ProteinFragment> result = new LinkedList<ProteinFragment>();
		result = collectFragments(query, clusters, extent);

		ProteinFragment resultFragment = assembleProtein(result, extent);
		return resultFragment;
	}
}
