package bioinfo.proteins.fragm3nt;

import java.util.ArrayList;

public class ClusterAnalysis {
	private ArrayList<FragmentCluster> clusters;
	
	public ClusterAnalysis(ArrayList<FragmentCluster> clusters) {
		this.clusters = new ArrayList<FragmentCluster>();
		for(FragmentCluster f: clusters) {
			this.clusters.add(f);
		}
	}

}
