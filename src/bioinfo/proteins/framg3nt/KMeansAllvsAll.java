package bioinfo.proteins.framg3nt;

import java.util.LinkedList;

import bioinfo.superpos.Kabsch;
import bioinfo.superpos.Transformation;

public class KMeansAllvsAll extends ClusterAlgorithm {

	public KMeansAllvsAll(LinkedList<ProteinFragment> f) {
		super(f);
	}

	@Override
	public void initializeClusters() {
		double[][][] kabschFood = new double[2][fragments.peek().fragLength][3];
		double minRMSD = Double.MAX_VALUE;
		Transformation t;
		ProteinFragment tmpFrag;
		ProteinFragment[] minPair = new ProteinFragment[2];

		// for every fragment pair...
		for (int i = 0; i < fragments.size(); i++) {
			if(fragments.get(i).getClusterIndex() > -1)
				continue;
			minRMSD = Double.MAX_VALUE;
			for (int j = 0; j < fragments.size(); j++) {
				// except of course every fragment with itself, that would
				// explode the experiment
				if (j == i)
					continue;
				// calculate the RMSD...
				kabschFood[0] = fragments.get(i).getAllResidues();
				kabschFood[1] = fragments.get(j).getAllResidues();
				t = Kabsch.calculateTransformation(kabschFood);
				
				// and find the pair with the minimal RMSD.
				double temp = t.getRmsd();
//				if(temp > 0)
//					System.err.println("here boss");
				if (minRMSD > temp) {
					minRMSD = temp;
					minPair = new ProteinFragment[2];
					minPair[0] = fragments.get(i);
					minPair[1] = fragments.get(j);
				}
			}

			// if one of the two fragments already belongs to a cluster, add the
			// other fragment to the same cluster
			if (minPair[1].getClusterIndex() > -1) {
				minPair[0].setClusterIndex(minPair[1].getClusterIndex());
				clusters.get(minPair[1].getClusterIndex()).add(minPair[0]);
			}
			// else create a new cluster and shove both fragments in this
			// cluster.
			else {
				clusters.addLast(new FragmentCluster());
				clusters.getLast().setCentroid(minPair[0]);
				minPair[0].setClusterIndex(clusters.size() - 1);
				minPair[1].setClusterIndex(clusters.size() - 1);
				clusters.getLast().add(minPair[0]);
				clusters.getLast().add(minPair[1]);
			}
		}
		System.out.println();
	}
}
