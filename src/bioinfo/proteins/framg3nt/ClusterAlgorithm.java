package bioinfo.proteins.framg3nt;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.util.LinkedList;

import org.zkoss.jsp.zul.FCKeditorTag;

import bioinfo.superpos.Kabsch;
import bioinfo.superpos.Transformation;

public abstract class ClusterAlgorithm {

	LinkedList<ProteinFragment> fragments = new LinkedList<ProteinFragment>();
	LinkedList<FragmentCluster> clusters = new LinkedList<FragmentCluster>();

	public ClusterAlgorithm(LinkedList<ProteinFragment> fr) {
		fragments = fr;
	}

	/**
	 * this method goes through the whole input and assigns or creates clusters.
	 * Always override.
	 */
	public void initializeClusters() {

	}

	/**
	 * this method has two components; it first calculates all cluster centroids
	 * and then reassigns all values to a cluster.
	 * 
	 * @return true if at least one updated cluster was changed (new centroid),
	 *         false if not.
	 */
	public boolean updateClusters() {
		boolean updated = false;
		updated = calculateCentroids();
		flushClusters();
		assignInstances();
		return updated;
	}

	private void flushClusters() {
		for (FragmentCluster c : clusters) {
			c.flush();
		}
	}

	public void assignInstances() {
		double[][][] kabschFood = new double[2][fragments.peek().fragLength][3];
		double minRMSD = Double.MAX_VALUE;
		FragmentCluster tempCluster = new FragmentCluster();
		Transformation t;
		int cur = -1;

		for (ProteinFragment f : fragments) {
			for (FragmentCluster cluster : clusters) {
				kabschFood[0] = f.getAllResidues();
				kabschFood[1] = cluster.getCentroid().getAllResidues();
				t = Kabsch.calculateTransformation(kabschFood);

				// and find the pair with the minimal RMSD.
				double temp = t.getRmsd();
				if (minRMSD > temp) {
					minRMSD = temp;
					tempCluster = cluster;
				}
			}
			cur = clusters.indexOf(tempCluster);
			clusters.get(cur).add(f);
		}

		for (FragmentCluster c : clusters) {
			System.out.println("============================");
			System.out.println(c.getCentroid().toString());
		}
		cur++;
	}

	/**
	 * calculates the new centroids for the cluster and compares them to the old
	 * ones
	 * 
	 * @return false if the clusters are the same (so the clusters are the same)
	 *         true if we need to reiterate
	 */
	public boolean calculateCentroids() {
		ProteinFragment curCentroid;
		boolean updated = false;
		for (FragmentCluster c : clusters) {
			curCentroid = new ProteinFragment(c.getCentroid().getID(), c
					.getCentroid().getAllResidues(), c.getCentroid()
					.getStartIndex(), c.getCentroid().getFragmentLength());
			c.calculateCentroid();
			if (!updated && c.equals(curCentroid))
				updated = true;
		}
		return updated;
	}

	/**
	 * this function performs the update operation on the clusters for n times
	 * or until the clusters aren't updated anymore - whatever happens first.
	 * 
	 * @param n
	 *            how many loops to run over the update function
	 */
	public void update(int n) {
		boolean updated = false;
		for (int i = 0; i < n; i++) {
			updated = updateClusters();
			if (!updated)
				break;
		}
	}

	/**
	 * this function performs the updateClusters operation until the clusters
	 * aren't updated anymore. Might be computationally expensive, if there are
	 * a lot of clusters
	 */
	public void update() {
		boolean updated = true;
		while (updated) {
			updated = updateClusters();
		}
	}

	public void toTextFiles(String prefix) {
		for (FragmentCluster c : clusters) {
			try {
				BufferedWriter br = new BufferedWriter(new FileWriter(prefix
						+ "_" + c.getCentroid().getID()));
				br.write(c.toString());
				br.close();
			} catch (Exception e) {
				e.printStackTrace();
			}
		}
	}
}
