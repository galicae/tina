package test;

import java.util.LinkedList;
import java.util.List;

import bioinfo.proteins.fragm3nt.ClusterAnalysis;
import bioinfo.proteins.fragm3nt.FragmentCluster;

public class AnalysisTest {

	public static void main(String[] args) {
		ClusterAnalysis ccl = new ClusterAnalysis("./clusters/");
		List<FragmentCluster> clusters = new LinkedList<FragmentCluster>();
		
		clusters = ccl.getClusters();
		double tempRMSD = 0;
		for(FragmentCluster c: clusters) {
			tempRMSD += ccl.calcClusterRMSD(c);
		}
		
		System.out.println("mean RMSD of kmeans: " + tempRMSD / clusters.size() + " in " + clusters.size() + " clusters");
	}

}
