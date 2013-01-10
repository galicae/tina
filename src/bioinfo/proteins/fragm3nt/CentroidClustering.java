package bioinfo.proteins.fragm3nt;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.util.LinkedList;

import bioinfo.superpos.Kabsch;
import bioinfo.superpos.Transformation;

/**
 * a class for centroid-based clustering methods, such as k-means and its
 * variations
 * 
 * @author nikos
 * 
 */
public abstract class CentroidClustering {
	LinkedList<ProteinFragment> fragments = new LinkedList<ProteinFragment>();
	LinkedList<FragmentCluster> clusters = new LinkedList<FragmentCluster>();

	public CentroidClustering(LinkedList<ProteinFragment> fr) {
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
		// System.out.println("starting update...");
		boolean updated = false;
		updated = calculateCentroids();
		flushClusters();
		assignInstances();
		// System.out.println("update closed");
		return updated;
	}

	/**
	 * empties clusters from their protein fragments
	 */
	private void flushClusters() {
		for (FragmentCluster c : clusters) {
			c.flush();
		}
	}

	/**
	 * goes over all protein fragments and assigns them to one of the clusters
	 * calculated by the init method
	 */
	public void assignInstances() {
		double[][][] kabschFood = new double[2][fragments.peek().fragLength][3];
		double minRMSD = Double.MAX_VALUE;
		FragmentCluster tempCluster = new FragmentCluster();
		Transformation t;
		// System.out.println("===================================\nCalling update\n===================================");
		// System.out.println(clusters.size());

		for (ProteinFragment f : fragments) {
			kabschFood[0] = f.getAllResidues();
			minRMSD = Double.MAX_VALUE;
			for (FragmentCluster cluster : clusters) {
				kabschFood[1] = cluster.getCentroid().getAllResidues();
				if (kabschFood[1].length == 5) {
					// System.out.println(f.getID() + " " +
					// cluster.getCentroid().getID());
					kabschFood[1] = cluster.getCentroid().getAllResidues();
				}
				t = Kabsch.calculateTransformation(kabschFood);

				// and find the pair with the minimal RMSD.
				double temp = t.getRmsd();
				if (minRMSD >= temp) {
					minRMSD = temp;
					tempCluster = cluster;
				}
			}
			// System.out.println("assign fragment " + f.getID() +
			// " to cluster " + tempCluster.getCentroid().getID() +
			// " with RMSD " + minRMSD);
			if (minRMSD < 2) {
				kabschFood[0] = tempCluster.getCentroid().getAllResidues();
				kabschFood[1] = f.getAllResidues();
				t = Kabsch.calculateTransformation(kabschFood);
				f.setCoordinates(t.transform(f.getAllResidues()));
				tempCluster.add(f);
			} else {
				clusters.addLast(new FragmentCluster());
				clusters.getLast().setCentroid(f);
				clusters.getLast().add(f);
			}
		}
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
		boolean needMoarUpdate = true;
		for (FragmentCluster c : clusters) {
			curCentroid = new ProteinFragment(c.getCentroid().getID(), c
					.getCentroid().getAllResidues(), c.getCentroid()
					.getStartIndex(), c.getCentroid().getFragmentLength());
			c.calculateCentroid();
			if (c.getCentroid().equals(curCentroid)) {
				needMoarUpdate = false;
				System.out.println("cluster " + c.getCentroid().getID()
						+ " has the same centroid!");
			} else
				needMoarUpdate = true;
		}
		return needMoarUpdate;
	}

	/**
	 * this function performs the update operation on the clusters for n times
	 * or until the clusters aren't updated anymore - whatever happens first.
	 * 
	 * @param n
	 *            how many loops to run over the update function
	 */
	public void update(int n) {
		System.out.println("Starting update....");
		boolean updated = true;
		for (int i = 0; i < n; i++) {
			System.out.println("iteration " + i);
			updated = updateClusters();
			for (FragmentCluster f : (LinkedList<FragmentCluster>) clusters
					.clone()) {
				if (f.getSize() == 0)
					clusters.remove(f);
			}
		}

	}

	/**
	 * this function performs the updateClusters operation until the clusters
	 * aren't updated anymore. Might be computationally expensive, if there are
	 * a lot of clusters
	 */
	public void update() {
		System.out.println("Starting update....");
		boolean updated = true;
		while (updated) {
			updated = updateClusters();
			for (FragmentCluster f : (LinkedList<FragmentCluster>) clusters
					.clone()) {
				if (f.getSize() == 0)
					clusters.remove(f);
			}
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

	public LinkedList<FragmentCluster> getClusters() {
		return this.clusters;
	}
}
