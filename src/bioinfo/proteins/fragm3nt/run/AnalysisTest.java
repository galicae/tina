package bioinfo.proteins.fragm3nt.run;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.util.LinkedList;
import java.util.List;

import bioinfo.proteins.fragm3nt.ClusterAnalysis;
import bioinfo.proteins.fragm3nt.FragmentCluster;

/**
 * cluster analysis class. Pretty self-explanatory.
 * 
 * @author galicae
 * 
 */
public class AnalysisTest {

	public static void main(String[] args) throws Exception {
		ClusterAnalysis ccl = new ClusterAnalysis("./clusters2/");
		List<FragmentCluster> clusters = new LinkedList<FragmentCluster>();
		BufferedWriter w1 = new BufferedWriter(new FileWriter(
				"/home/galicae/Documents/intraCluster2"));

		clusters = ccl.getClusters();
		double tempRMSD = 0;
		double cur = 0;
		for (FragmentCluster c : clusters) {
			cur = ccl.calcIntraClusterRMSD(c);
			w1.write(cur + "\n");
			tempRMSD += cur;
		}
		w1.close();
		double intra = tempRMSD / clusters.size();
		double inter = ccl.calcInterClusterRMSD();

		System.out.println("clusters: " + clusters.size());
		System.out.println("mean intracluster RMSD: " + intra);
		System.out.println("mean intercluster RMSD " + inter);
	}

}
