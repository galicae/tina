package bioinfo.proteins.fragm3nt.assembler;

import java.util.LinkedList;

import bioinfo.superpos.*;

import bioinfo.alignment.matrices.QuasarMatrix;
import bioinfo.proteins.PDBEntry;
import bioinfo.proteins.fragm3nt.FragmentCluster;
import bioinfo.proteins.fragm3nt.ProteinFragment;

/**
 * a cheat version of Assembler. The fragment selection method is overridden:
 * instead of choosing based on sequence similarity (original Assembler),
 * CheatAssembler registers the maximum sequence score of a fragment (as well as
 * its minimum) and chooses the fragment that best superposes on the native
 * structure. Thus the assembler allows to evaluate whether the PSSM method is
 * adequate as a fragment selection method (hint: it is not; it is woefully
 * bad).
 * 
 * @author galicae
 * 
 */
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

	/**
	 * giant function - extremely ugly, but time was pressing. This function
	 * wasn't meant to be used as a standard for the server, and it shouldn't.
	 * Created for a "quick test" that quickly escalated to full-blown
	 * evaluation work
	 * 
	 * @param clusters
	 *            the cluster library
	 * @param extent
	 *            to which extent to align fragments with one another (overlap)
	 * @param query
	 *            the sequence of the protein in question
	 * @return the assembled prediction
	 */
	@Deprecated
	private ProteinFragment assembleProtein(
			LinkedList<FragmentCluster> clusters, int extent, String query) {

		ProteinFragment tempResult = new ProteinFragment("bla",
				new double[0][0], extent);

		// get the native protein in ProteinFragment form, it is more
		// flexible this way
		ProteinFragment controlFrag = new ProteinFragment("control",
				PDBReduce.reduceSinglePDB(control), extent);
		controlFrag.setSequence(query);

		ProteinFragment resultFragment = new ProteinFragment("bla",
				new double[0][0], extent);
		// for the first fragment cheat and just copy it from the native
		resultFragment = controlFrag.getPart(0, fragLength);
		int size = clusters.size();

		// after the first fragment is snatched up we have predicted till index
		// fragmentLength
		int index = fragLength;
		// define the variables we need: we want to monitor whether the sequence
		// score of the fragments we pick up actually coincides with the maximum
		// sequence score available, that is with the fragment that the method
		// would normally pick up.
		int indexLimit = query.length() - fragLength + extent;
		double minClusterRMScore = Double.MAX_VALUE; // we want the cluster
		// that provides the minimum RMSD from the native structure
		double maxClusterSeqScore = Integer.MIN_VALUE; // we also want to
		// capture the fragment that would be selected by the "normal" method

		double minClusterSeqScore = Integer.MAX_VALUE; // we also want to
		// capture the fragment that would be selected by the "normal" method
		double curClusterSeqScore = 0; // temporary variable for the sequence
										// score
		double actualClusterSeqScore = 0; // this stores the value of the
		// sequence score of the fragment that minimizes RMSD
		double actualScore = 0; // the cumulative score of the fragments picked
								// up
		double maxScore = 0; // the cumulative score of the fragments the
		// method would normally pick up
		double minScore = 0;

		double curClusterRMScore = 0; // temporary variable for the RMSD score
		// between a cluster and the predicted prefix

		double[][] matrix = QuasarMatrix.DAYHOFF_MATRIX;

		// so, while we're not at the end (index is how much we've already
		// predicted)...
		while (index < indexLimit) {
			ProteinFragment tempResultFragment = clusters.get(0).getCentroid()
					.clone();

			minClusterRMScore = Double.MAX_VALUE;
			maxClusterSeqScore = Integer.MIN_VALUE;
			minClusterSeqScore = Integer.MAX_VALUE;

			// ...go through all clusters...
			for (int i = 0; i < size; i++) {

				// ...calculate the sequence score of the current pairing...
				curClusterSeqScore = 0;

				FragmentCluster c = new FragmentCluster();
				c = clusters.get(i);
				for (int j = 0; j < fragLength; j++) {
					int y = query.charAt(j) - 65;
					for (int k = 0; k < 25; k++) {
						curClusterSeqScore += matrix[y][k] * c.getPssm()[j][k];
					}
				}

				// always clone stuff because you are dealing with a mixture
				// of primitive data types and wrappers
				tempResultFragment = resultFragment.clone();
				c.getCentroid().setSequence(
						query.substring(index - extent, index - extent
								+ fragLength));
				alignFragments(tempResultFragment, clusters.get(i)
						.getCentroid(), extent);
				// ...calculate the RMSD of the current pair...
				curClusterRMScore = compare(
						controlFrag.getPart(0,
								tempResultFragment.getAllResidues().length),
						tempResultFragment);

				// save the maximum sequence score found in the round
				if (maxClusterSeqScore < curClusterSeqScore) {
					maxClusterSeqScore = curClusterSeqScore;
				}
				if (minClusterSeqScore > curClusterSeqScore) {
					minClusterSeqScore = curClusterSeqScore;
				}

				// save minimum RMSD of the round and the sequence score
				// of the same fragment
				if (minClusterRMScore > curClusterRMScore) {
					minClusterRMScore = curClusterRMScore;
					tempResult = tempResultFragment.clone();
					actualClusterSeqScore = curClusterSeqScore;
				}
			}

			// increment the cumulative scores
			actualScore += actualClusterSeqScore;
			maxScore += maxClusterSeqScore;
			minScore += minClusterSeqScore;
			// and save the best-scoring fragment
			resultFragment = tempResult.clone();
			// update index or you get an infinite loop
			index = resultFragment.getAllResidues().length;
		}

		// repeat for the last fragment
		int newExtent = fragLength
				- (query.length() - resultFragment.getAllResidues().length);
		for (int i = 0; i < size; i++) {

			curClusterSeqScore = 0;
			FragmentCluster c = new FragmentCluster();
			c = clusters.get(i);
			for (int j = 0; j < fragLength; j++) {
				int y = query.charAt(j) - 65;
				for (int k = 0; k < 25; k++) {
					curClusterSeqScore += matrix[y][k] * c.getPssm()[j][k];
				}
			}

			ProteinFragment temp = clusters.get(0).getCentroid().clone();
			temp = resultFragment.clone();
			clusters.get(i)
					.getCentroid()
					.setSequence(
							query.substring(index - newExtent, index
									- newExtent + fragLength));
			alignFragments(temp, clusters.get(i).getCentroid(), newExtent);
			curClusterRMScore = compare(
					controlFrag.getPart(0, temp.getAllResidues().length), temp);
			if (maxClusterSeqScore < curClusterSeqScore) {
				maxClusterSeqScore = curClusterSeqScore;
			}
			if (minClusterSeqScore < curClusterSeqScore) {
				minClusterSeqScore = curClusterSeqScore;
			}
			if (minClusterRMScore > curClusterRMScore) {
				minClusterRMScore = curClusterRMScore;
				tempResult = temp.clone();
				actualClusterSeqScore = curClusterSeqScore;
			}
		}
		actualScore += actualClusterSeqScore;
		maxScore += maxClusterSeqScore;
		minScore += minClusterSeqScore;
		resultFragment = tempResult.clone();
		index = resultFragment.getAllResidues().length;

		// print results. This is not really clean, but I see no other way
		// that wouldn't involve major restructuring of ProteinFragment
		// and there's no need to do that, the class is already overblown
		System.out
				.print(actualScore + "\t" + maxScore + "\t" + minScore + "\t");
		resultFragment.setSequence(query);
		return resultFragment;
	}

	/**
	 * this function calculates the RMSD between two ProteinFragments of the
	 * same size
	 * 
	 * @param nativStr
	 *            the first ProteinFragment
	 * @param predStr
	 *            the second ProteinFragment
	 * @return the RMSD of the optimal superposition of fragment predStr to
	 *         nativeStr
	 */
	private double compare(ProteinFragment nativStr, ProteinFragment predStr) {
		double[][][] kabschFood = new double[2][predStr.getAllResidues().length][3];
		kabschFood[0] = nativStr.getAllResidues();
		kabschFood[1] = predStr.getAllResidues();
		Transformation t = Kabsch.calculateTransformation(kabschFood);
		return t.getRmsd();
	}

}
