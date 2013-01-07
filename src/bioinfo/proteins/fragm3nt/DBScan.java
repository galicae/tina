package bioinfo.proteins.fragm3nt;

import java.util.ArrayList;
import java.util.LinkedList;

import bioinfo.superpos.Kabsch;
import bioinfo.superpos.Transformation;

/**
 * an implementation of the DBSCAN (density-based clustering) algorithm
 * described in
 * "A density-based algorithm for discovering clusters in large spatial databases with noise"
 * by Martin Ester, Hans-Peter Kriegel, J�rg Sander, Xiaowei Xu (1996-)
 * 
 * mostly based on the wikipedia pseudocode of the algorithm
 * 
 * @author nikos
 * 
 */
public class DBScan {
	int MINPTS;
	double EPS;
	// the distance matrix needn't be quadratic: we assume that dist(x,y) =
	// dist(y,x), and since dist(x,x)=0 we only need to save n(n-1)/2 distances.
	// Thus the matrix must only be (almost) half as big. We need a smart index
	// structure to access the information that would be stored in (x,y).
	int[] distance;

	// since we use a triangular matrix it is only natural that we use the
	// triangular numbers. It did take some time until I was aided to this
	// discovery, but since j is always larger than i we can reach point (i,j)
	// by going to j(j-1)/2 + i

	/**
	 * this function computes the pairwise distances of all protein fragments in
	 * a list and saves them in an one-dimensional array, making the assumptions
	 * listed above
	 * 
	 * @param data
	 *            the list of protein fragments
	 */
	private void calculateAllDistances(ArrayList<ProteinFragment> data) {
		double[][][] kabschFood = new double[2][data.get(0).fragLength][3];
		Transformation t;
		int temp = 0;
		distance = new int[((data.size() - 1) * data.size()) / 2];
		for (int i = 0; i < data.size() - 1; i++) {
			for (int j = i + 1; j < data.size(); j++) {
				// don't allow a fragment to be compared with itself. Computing
				// time is precious!
				if (j == i)
					continue;
				kabschFood[0] = data.get(i).getAllResidues();
				kabschFood[1] = data.get(j).getAllResidues();
				t = Kabsch.calculateTransformation(kabschFood);
				temp = (int) (t.getRmsd() * 1000);
				distance[((j * (j - 1)) / 2) + i] = temp;
			}
		}
	}


	/**
	 * finds all neighbours of fragment at position 'position', defined by their
	 * distance to said fragment. If there are more neighbours than the
	 * predefined EPS cut-off, the neighbours are returned.
	 * 
	 * @param position
	 *            the index of the protein fragment
	 * @param data
	 *            the list of protein fragments
	 * @return null, if there are not enough neighbours, else the neighbours
	 *         themselves.
	 */
	private ArrayList<ProteinFragment> getNeighbours(int position,
			ArrayList<ProteinFragment> data) {
		ArrayList<ProteinFragment> neighbPoints = new ArrayList<ProteinFragment>();
		for (int i = 0; i < data.size(); i++) {
			if (position == i)
				continue;
			if (getDistAt(position, i) < EPS) {
				neighbPoints.add(data.get(i));
			}
		}

		return neighbPoints;
	}

	/**
	 * shortcut to retrieve a two-dimensional position from the one-dimensional
	 * array.
	 * 
	 * @param i
	 *            the x-coordinate
	 * @param j
	 *            the y-coordinate
	 * @return
	 */
	private int getDistAt(int i, int j) {
		if (j > i)
			return distance[((j * (j - 1)) / 2) + i];
		else
			return distance[((i * (i - 1)) / 2) + j];
	}

	/**
	 * this function expands clusters by checking the neighbours of a fragment
	 * for the density condition
	 * 
	 * @param p
	 *            the protein fragment whose neighbourhood we are checking
	 * @param c
	 *            the cluster to which p belongs
	 * @param neighbours
	 *            the direct neighbours of p
	 * @param data
	 *            the whole fragment library
	 */
	public void expandCluster(ProteinFragment p, FragmentCluster c,
			ArrayList<ProteinFragment> neighbours,
			ArrayList<ProteinFragment> data) {
		c.add(p);
		for (ProteinFragment f: neighbours) {			
			if (!f.isVisited()) {
				f.setVisited(true);
				ArrayList<ProteinFragment> fNeighbours = getNeighbours(
						data.indexOf(f), data);
				if (fNeighbours.size() >= MINPTS)
					neighbours.addAll(fNeighbours);
			}
			if (f.getClusterIndex() == -1) {
				f.setClusterIndex(c.getCentroid().getClusterIndex());
				c.add(f);
			}
		}
	}
	
	
	public void oppaDBStyle(int minpts, double eps, ArrayList<ProteinFragment> data, LinkedList<FragmentCluster> clusters) {
		MINPTS = minpts;
		EPS = eps;
		ProteinFragment p = new ProteinFragment(null, null, 0, 0);
		for(int i = 0; i < data.size(); i++) {
			p = data.get(i);
			if(!p.isVisited()) {
				p.setVisited(true);
				ArrayList<ProteinFragment> pNeighbours = getNeighbours(i, data);
				if(pNeighbours.size() < MINPTS) {
					p.setNoise(true);
				}
				else {
					clusters.addLast(new FragmentCluster());
					p.setClusterIndex(clusters.size() - 1);
					clusters.getLast().setCentroid(p);
					expandCluster(p, clusters.getLast(), pNeighbours, data);
				}
			}
		}
	}
}
