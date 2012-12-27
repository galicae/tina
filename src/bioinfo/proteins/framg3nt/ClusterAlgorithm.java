package bioinfo.proteins.framg3nt;

import java.util.LinkedList;

public abstract class ClusterAlgorithm {
	
	LinkedList<ProteinFragment> fragments = new LinkedList<ProteinFragment>();
	LinkedList<FragmentCluster> clusters = new LinkedList<FragmentCluster>();
	
	

	public ClusterAlgorithm(LinkedList<ProteinFragment> fr) {
		fragments = fr;
	}
	
	
	/**
	 * this method goes through the whole input and assigns or creates clusters. Always override.
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
		assignInstances();
		return updated;
	}

	public void assignInstances() {

	}

	/**
	 * calculates the new centroids for the cluster and compares them to the old ones
	 * @return false if the clusters are the same (so the clusters are the same) true if we need to reiterate
	 */
	public boolean calculateCentroids() {
		ProteinFragment curCentroid;
		boolean updated = false;
		for(FragmentCluster c: clusters) {
			curCentroid = c.getCentroid();
			c.calculateCentroid();
			if(!updated && c.equals(curCentroid))
				updated = true;
		}
		return updated;
	}
}
