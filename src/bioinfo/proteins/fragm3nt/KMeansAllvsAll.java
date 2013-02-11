package bioinfo.proteins.fragm3nt;

import java.util.LinkedList;

import bioinfo.superpos.Kabsch;
import bioinfo.superpos.Transformation;

public class KMeansAllvsAll extends CentroidClustering {

	public KMeansAllvsAll(LinkedList<ProteinFragment> f, double accuracy) {
		super(f, accuracy);
	}

	/**
	 * this function initializes clusters to use the k-means clustering
	 * algorithm. The distance function here is the RMSD of the fragments to one
	 * another, or to a cluster centroid. This initialization function is a
	 * slight modification of the classical k-means algorithm, in that it
	 * doesn't assign a fixed number of clusters but rather puts each fragment
	 * in a cluster with its closest neighbor. This is a definite overestimation
	 * of the number of clusters, but nevertheless beats the assignment of k,
	 * where k=random.
	 */
	@Override
	public void initializeClusters(int p) {
		if (p > 0) {
			double[][][] kabschFood = new double[2][fragments.peek().fragLength][3];
			double minRMSD = Double.MAX_VALUE;
			Transformation t;
			ProteinFragment[] minPair = new ProteinFragment[2];

			// error filtering
			LinkedList<ProteinFragment> correctList = new LinkedList<ProteinFragment>();
			LinkedList<ProteinFragment> wrongList = new LinkedList<ProteinFragment>();
			for (ProteinFragment f : fragments) {
				if (checkFragment(f))
					correctList.add(f);
				else
					wrongList.add(f);
			}
			System.out.println("Sorted fragments in good and bad");
			fragments = correctList;

			for (ProteinFragment f : wrongList) {
				f.getSequence();
			}
			// for every fragment pair...
			for (int i = 0; i < fragments.size(); i++) {
				if (fragments.get(i).getClusterIndex() > -1)
					continue;
				minRMSD = Double.MAX_VALUE;
				kabschFood[0] = fragments.get(i).getAllResidues();
				for (int j = 0; j < fragments.size(); j++) {
					// except of course every fragment with itself, that would
					// explode the experiment
					if (j == i)
						continue;
					// calculate the RMSD...

					kabschFood[1] = fragments.get(j).getAllResidues();
					t = Kabsch.calculateTransformation(kabschFood);

					// and find the pair with the minimal RMSD.
					double temp = t.getRmsd();
					// if(temp > 0)
					// System.err.println("here boss");
					if (minRMSD > temp) {
						minRMSD = temp;
						minPair = new ProteinFragment[2];
						minPair[0] = fragments.get(i);
						minPair[1] = fragments.get(j);
					}
				}

				// if one of the two fragments already belongs to a cluster, add
				// the
				// other fragment to the same cluster
				if (minPair[1].getClusterIndex() > -1) {
					minPair[0] = minPair[0].clone();
					minPair[1] = minPair[1].clone();
					minPair[0].setClusterIndex(minPair[1].getClusterIndex());

					kabschFood[0] = clusters.get(minPair[0].getClusterIndex())
							.getCentroid().getAllResidues();
					kabschFood[1] = minPair[0].getAllResidues();
					t = Kabsch.calculateTransformation(kabschFood);
					minPair[0].setCoordinates(t.transform(minPair[0]
							.getAllResidues()));
					clusters.get(minPair[1].getClusterIndex()).add(minPair[0]);
				}
				// else create a new cluster and shove both fragments in this
				// cluster.
				else {
					minPair[0] = minPair[0].clone();
					minPair[1] = minPair[1].clone();
					clusters.addLast(new FragmentCluster());
					clusters.getLast().setCentroid(minPair[1]);
					minPair[0].setClusterIndex(clusters.size() - 1);
					minPair[1].setClusterIndex(clusters.size() - 1);
					kabschFood[0] = minPair[1].getAllResidues();
					kabschFood[1] = minPair[0].getAllResidues();
					t = Kabsch.calculateTransformation(kabschFood);
					minPair[0].setCoordinates(t.transform(minPair[0]
							.getAllResidues()));
					clusters.getLast().add(minPair[0]);
					clusters.getLast().add(minPair[1]);
				}
			}
		} else {
			for (ProteinFragment f : fragments) {
				clusters.addLast(new FragmentCluster());
				clusters.getLast().setCentroid(f);
				clusters.getLast().add(f);
			}
		}
	}

}
