package bioinfo.proteins.fragm3nt;

import java.util.ArrayList;

import bioinfo.superpos.Kabsch;
import bioinfo.superpos.Transformation;

/**
 * an implementation of the DBSCAN (density-based clustering) algorithm
 * described in
 * "A density-based algorithm for discovering clusters in large spatial databases with noise"
 * by Martin Ester, Hans-Peter Kriegel, Jörg Sander, Xiaowei Xu (1996-)
 * 
 * mostly based on the wikipedia pseudocode of the algorithm
 * 
 * @author nikos
 * 
 */
public class DBScan {
	int MINPTS;
	double EPS;
	ArrayList<ProteinFragment> DATA = new ArrayList<ProteinFragment>();
	// the distance matrix needn't be quadratic: we assume that dist(x,y) =
	// dist(y,x), and since dist(x,x)=0 we only need to save n(n-1)/2 distances.
	// Thus the matrix must only be (almost) half as big. We need a smart index
	// structure to access the information that would be stored in (x,y).
	int[] distance;

	// since we use a triangular matrix it is only natural that we use the
	// triangular numbers. It did take some time until I was aided to this
	// discovery, but since j is always larger than i we can reach point (i,j)
	// by going to j(j-1)/2 + i

	public void calculateAllDistances() {
		double[][][] kabschFood = new double[2][DATA.get(0).fragLength][3];
		Transformation t;
		distance = new int[((DATA.size() - 1) * DATA.size()) / 2];
		for (int i = 0; i < DATA.size() - 1; i++) {
			for (int j = i + 1; j < DATA.size(); j++) {
				// don't allow a fragment to be compared with itself. Computing
				// time is precious!
				if (j == i)
					continue;
				kabschFood[0] = DATA.get(i).getAllResidues();
				kabschFood[1] = DATA.get(j).getAllResidues();
				t = Kabsch.calculateTransformation(kabschFood);
				distance[((j * (j - 1)) / 2) + i] = (int) (t.getRmsd() * 1000);
			}
		}
	}
	
//	DBSCAN(D, eps, MinPts)
//	   C = 0
//	   for each unvisited point P in dataset D
//	      mark P as visited
//	      NeighborPts = regionQuery(P, eps)
//	      if sizeof(NeighborPts) < MinPts
//	         mark P as NOISE
//	      else
//	         C = next cluster
//	         expandCluster(P, NeighborPts, C, eps, MinPts)
//	          
//	expandCluster(P, NeighborPts, C, eps, MinPts)
//	   add P to cluster C
//	   for each point P' in NeighborPts 
//	      if P' is not visited
//	         mark P' as visited
//	         NeighborPts' = regionQuery(P', eps)
//	         if sizeof(NeighborPts') >= MinPts
//	            NeighborPts = NeighborPts joined with NeighborPts'
//	      if P' is not yet member of any cluster
//	         add P' to cluster C
//	          
//	regionQuery(P, eps)
//	   return all points within P's eps-neighborhood

}
