package bioinfo.proteins.fragm3nt;

import java.util.ArrayList;
import java.util.LinkedList;
import java.util.List;

import bioinfo.proteins.PDBEntry;
import bioinfo.proteins.PDBFileReader;

public class CollectiveClustering {
	private int fragLength;
	private String proteins;
	
	public CollectiveClustering(int fragLength, String proteins) {
		this.fragLength = fragLength;
		this.proteins = proteins;
	}
	
	/**
	 * collective function for kmeans clustering. Seems to work.
	 * @param update
	 * @return
	 */
	public LinkedList<FragmentCluster> runKmeans(int update) {
		PDBFileReader reader = new PDBFileReader(proteins);
		List<PDBEntry> files = new LinkedList<PDBEntry>();
		LinkedList<ProteinFragment> pList = new LinkedList<ProteinFragment>();
		files = reader.readPdbFolder();
		// files.add(pdb1);
		for (PDBEntry e : files) {
			Fragmenter.crunchBackboneSeq(e, pList, fragLength);
		}
		KMeansAllvsAll clustah = new KMeansAllvsAll(pList);
		LinkedList<FragmentCluster> clusters = new LinkedList<FragmentCluster>();
		clustah.initializeClusters();
		clustah.update(update);

		clusters = clustah.getClusters();
		return clusters;
	}
	
	public LinkedList<FragmentCluster> runDBScan(int minpts, double eps) {
		PDBFileReader reader = new PDBFileReader(proteins);
		// and save it in a files list
		List<PDBEntry> files = new LinkedList<PDBEntry>();
		ArrayList<ProteinFragment> pList = new ArrayList<ProteinFragment>();
		
		files = reader.readPdbFolder();
		
		// next make fragments out of all PDBs and save them in pList
		for(PDBEntry e: files) {
			Fragmenter.crunchBackboneSeq(e, pList, 5);
		}
		DBScan clustah = new DBScan();
		LinkedList<FragmentCluster> clusters = new LinkedList<FragmentCluster>();
		clustah.oppaDBStyle(minpts, eps, pList, clusters);
		return clusters;
	}
}
