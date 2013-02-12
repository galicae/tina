package test;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.util.LinkedList;
import java.util.List;

import bioinfo.proteins.PDBEntry;
import bioinfo.proteins.PDBFileReader;
import bioinfo.proteins.fragm3nt.FragmentCluster;
import bioinfo.proteins.fragm3nt.Fragmenter;
import bioinfo.proteins.fragm3nt.KMeansAllvsAll;
import bioinfo.proteins.fragm3nt.ProteinFragment;

public class KMeansMain {
	public static void main(String[] args) {
		int updCycles = Integer.parseInt(args[0]);
		int fragLength = Integer.parseInt(args[1]);
		double acc = Double.parseDouble(args[2]);
		String pre = args[3];

		long startTime = System.currentTimeMillis();
		long readTime = 0;
		long crunchTime = 0;
		long initTime = 0;
		long updateTime = 0;
		long writeTime = 0;
		long finishTime = 0;
		PDBFileReader reader = new PDBFileReader(
				"/home/p/papadopoulos/Desktop/proteins");
		List<PDBEntry> files = new LinkedList<PDBEntry>();

		System.out
				.println("=========================================================");
		System.out.println("KMeans clustering with fragment length of "
				+ fragLength + " and " + updCycles + " update cycles.");

		readTime = System.currentTimeMillis();
		System.out.println("reading ./proteins...");
		files = reader.readPdbFolder();
		crunchTime = System.currentTimeMillis();

		System.out.println("Reading performed in (ms) "
				+ (crunchTime - readTime));

		LinkedList<ProteinFragment> pList = new LinkedList<ProteinFragment>();
		System.out.println("crunching records...");
		for (PDBEntry e : files) {
			// TODO: use reeeeal function
			Fragmenter.crunchBackboneSeq(e, pList, fragLength);
		}
		initTime = System.currentTimeMillis();

		System.out.println("Crunched fragments in (ms) "
				+ (initTime - crunchTime));

		KMeansAllvsAll clusterAlgorithm = new KMeansAllvsAll(pList, acc);
		System.out.println("initializing clusters...");
		clusterAlgorithm.initializeClusters(1);
		// clusterAlgorithm.checkAllFragments();
		updateTime = System.currentTimeMillis();

		System.out.println("Cluster initialization in (ms) "
				+ (updateTime - initTime));
		System.out.println("updating " + clusterAlgorithm.getClusters().size()
				+ " clusters...");
		clusterAlgorithm.update(updCycles);
		// clusterAlgorithm.checkAllFragments();
		writeTime = System.currentTimeMillis();
		System.out.println("Updated " + updCycles + " times in (ms) "
				+ (writeTime - updateTime));

		// System.out.println(clusterAlgorithm.getClusters().getFirst().toString());

		for (FragmentCluster c : clusterAlgorithm.getClusters()) {
			try {
				BufferedWriter br = new BufferedWriter(new FileWriter("./"
						+ pre + "_" + c.getCentroid().getID()));
				br.write(c.toString());
				br.close();
			} catch (Exception e) {
				e.printStackTrace();
			}
		}
		finishTime = System.currentTimeMillis();
		System.out.println("Printed in (ms) " + (finishTime - writeTime));
		System.out
				.println("=========================================================");
		System.out.println("Started at " + startTime);
		System.out.println("Finished at " + finishTime);
		System.out.println("Total time of " + (finishTime - startTime)
				+ " (ms)");
	}
}
