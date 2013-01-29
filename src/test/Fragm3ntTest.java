package test;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.util.LinkedList;

import bioinfo.proteins.fragm3nt.CollectiveClustering;
import bioinfo.proteins.fragm3nt.FragmentCluster;

public class Fragm3ntTest {
	public static void main(String[] args) {
		LinkedList<FragmentCluster> clusters = new LinkedList<FragmentCluster>();
		CollectiveClustering db = new CollectiveClustering(5, "./proteins2");
		clusters = db.runKmeans(50, 0.5);
//		clusters = db.runDBScan(4, 1.0);
		
		int sumOfFrags = 0;
		for(FragmentCluster c: clusters) {
			sumOfFrags += c.getSize();
		}
		
		// report clustering efficiency
		System.out.format("%d fragments in %d clusters.\n" , sumOfFrags, clusters.size());
		
		// write out clusters
		for (FragmentCluster c : clusters) {
			try {
				BufferedWriter br = new BufferedWriter(new FileWriter("./clusters/k" + c.getCentroid().getID()));
				br.write(c.toString());
				br.close();
			} catch (Exception e) {
				e.printStackTrace();
			}
		}
	}
}

