package bioinfo.proteins.fr4gment;

import java.io.File;
import java.io.FileFilter;
import java.util.HashMap;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import bioinfo.proteins.PDBEntry;

/**
 * 
 * this class contains complete multiple superposition clustering
 * formatted Marie-style
 * peformed as clustring over pairwise superpositions and clustered over TM score in a SingleLinkage manner
 * "Kids, don't try this at home!" =)
 * @author andreseitz
 *
 */
public class MSClustering {
	
	private HashMap<String,Integer> pdb2cluster = new HashMap<String,Integer>();
	private HashMap<Integer, MultipleSuperposition> cluster2data = new HashMap<Integer, MultipleSuperposition>();
	private HashMap<Integer,Integer> cluster2fold = new HashMap<Integer, Integer>();
	private HashMap<Integer,int[]> fold2clusters = new HashMap<Integer, int[]>();
	private double threshold;
	
	/**
	 * reads complete clustering from Marie sorted in folder structure
	 * parameter is the folder containing all folds etc...
	 * folder strcuture has to look like:
	 * cluster05(:=foldername)
	 * 		|-1 (fold)
	 * 			|-20 (cluster)
	 * 			|-21 (cluster)
	 * 			...
	 * 		|-10 (fold)
	 * 			|...
	 * threshold describes the TM-score x thta was used to cut the pairwise superposition list with which the clustering was performed
	 * @param foldername
	 */
	public MSClustering(String foldername, double threshold){
		this.threshold = threshold;
		int foldnum;
		int[] clustersPerFold;
		int clusterCount;
		int clusterid;
		Pattern clusterNamePattern = Pattern.compile("cluster(\\d+)\\.pdb");
		Matcher clusterNameMatcher;
		File superfolder = new File(foldername);
		File[] folds = superfolder.listFiles(new FileFilter() {
			@Override
			public boolean accept(File pathname) {
				return pathname.isDirectory();
			}
		});
		File[] clusters;
		for(File fold : folds){
			foldnum = Integer.parseInt(fold.getName());
			clusters = fold.listFiles(new FileFilter() {
				@Override
				public boolean accept(File pathname) {
					return pathname.getName().endsWith(".pdb") ? true : false;
				}
			});
			clustersPerFold = new int[clusters.length];
			clusterCount = 0;
			for (File cluster : clusters){
				clusterNameMatcher = clusterNamePattern.matcher(cluster.getName());
				if(clusterNameMatcher.find()){
					clusterid = Integer.parseInt(clusterNameMatcher.group(1));
				}else{
					clusterid = -1;
					System.err.println(cluster.getName());
				}
				clustersPerFold[clusterCount] = clusterid;
				cluster2fold.put(clusterid, foldnum);
				cluster2data.put(clusterid, new MultipleSuperposition(cluster.getAbsolutePath()));
				for(PDBEntry entry : cluster2data.get(clusterid).getStructures()){
					pdb2cluster.put(entry.getID(),clusterid);
				}
				clusterCount++;
			}
			fold2clusters.put(foldnum, clustersPerFold);
		}
		System.out.println();
	}

	public static void main(String[] args) {
		MSClustering ms = new MSClustering("/Users/andreseitz/Desktop/voroeval/msp_tmp/cutoff05", 0);
	}
}
