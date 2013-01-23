package test;

import java.util.ArrayList;

import bioinfo.proteins.fragm3nt.ClusterAnalysis;
import bioinfo.proteins.fragm3nt.FragmentCluster;

public class AnalysisTest {

	public static void main(String[] args) {
		ClusterAnalysis ccl = new ClusterAnalysis("./clusters/");
		ArrayList<FragmentCluster> clusters = new ArrayList<FragmentCluster>();
		
		clusters = ccl.getClusters();
		System.out.println();
	}

}
