package bioinfo.proteins.fragm3nt;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.InputStreamReader;
import java.util.LinkedList;

import bioinfo.proteins.Atom;
import bioinfo.superpos.Kabsch;
import bioinfo.superpos.Transformation;

public class ClusterAnalysis {
	private LinkedList<FragmentCluster> clusters;

	public ClusterAnalysis(LinkedList<FragmentCluster> clusters) {
		this.clusters = new LinkedList<FragmentCluster>();
		for (FragmentCluster f : clusters) {
			this.clusters.add(f);
		}
	}

	/**
	 * this constructor takes care of the reading of the clusters.
	 * 
	 * @param folder
	 *            the path to the folder where the cluster files are saved
	 */
	public ClusterAnalysis(String folder) {
		this.clusters = new LinkedList<FragmentCluster>();
		BufferedReader br = null;
		try {
			if (folder == null) {
				System.err.println("null folder - no clusters read");
			}
			if (!folder.endsWith("/"))
				folder += "/";
			File[] files = new File(folder).listFiles();
			if (files.length == 0) {
				System.err.println("no files found - no clusters read");
			}
			for (File file : files) {
				br = new BufferedReader(new InputStreamReader(
						new FileInputStream(file)));
				clusters.add(parseCluster(br, file.getName()));
			}
		} catch (Exception e) {
			e.printStackTrace();
		}
	}

	// public ClusterAnalysis(LinkedList<FragmentCluster> clusters) {
	//
	// }

	/**
	 * this function takes care of the actual parsing of cluster files, as
	 * produced by the Fragm3nt pipeline
	 * 
	 * @param br
	 *            a buffered reader (should contain a FileReader)
	 * @param clusterId
	 *            the id of the cluster. Should be the name of the file without
	 *            any prefix
	 * @return the fragment cluster contained in the text file contained in br
	 */
	private FragmentCluster parseCluster(BufferedReader br, String clusterId) {
		try {
			String tempId = clusterId;
			FragmentCluster result = new FragmentCluster();
			result.setName(clusterId);
			String line = br.readLine();
			String seq = "";
			int fragLength = 0;
			LinkedList<Atom> atoms = new LinkedList<Atom>();
			String name;
			double[] coord = new double[3];
			while (line != null) {
				if (line.startsWith("REMARK")) {
					seq = line.split(" ")[2];
				} else if (line.startsWith("MODEL")) {
					tempId += "_" + line.split("\\s+")[1];
				} else if (line.startsWith("ATOM")) {
					name = line.substring(12, 16).trim();
					coord[0] = Double
							.parseDouble(line.substring(30, 38).trim());
					coord[1] = Double
							.parseDouble(line.substring(38, 46).trim());
					coord[2] = Double
							.parseDouble(line.substring(46, 54).trim());
					atoms.add(new Atom(name, coord));
					name = "";
					coord = new double[3];
				} else if (line.startsWith("ENDMDL")) {
					fragLength = atoms.size();
					Atom[] atomLinked = new Atom[atoms.size()];
					for (int i = 0; i < atoms.size(); i++) {
						Atom tempAt = atoms.get(i);
						atomLinked[i] = new Atom(tempAt.getType(),
								tempAt.getPosition());
					}
					ProteinFragment heyIJustMetYou = new ProteinFragment(
							tempId, seq, atomLinked, fragLength);
					// result.add(heyIJustMetYou);
					result.add(heyIJustMetYou);
					tempId = clusterId;
					atoms.clear();
				}

				line = br.readLine();
			}
			result.calculateCentroid();
			return result;
		} catch (Exception e) {
			e.printStackTrace();
		}

		return null;
	}

	public LinkedList<FragmentCluster> getClusters() {
		return this.clusters;
	}

	/**
	 * calculates the mean of the RMS distances of all fragments in a cluster
	 * 
	 * @param c
	 *            the cluster in question
	 * @return the mean of the RMS distances of all fragments
	 */
	public double calcIntraClusterRMSD(FragmentCluster c) {
		double[][][] kabschFood = new double[2][c.getFragmentLength()][3];
		Transformation t;
		double temp = 0;
		int size = c.getSize();
		// first calculate all distances
		for (int i = 0; i < size - 1; i++) {
			for (int j = i + 1; j < size; j++) {
				if (j == i)
					continue;
				kabschFood[0] = c.getFragments().get(i).getAllResidues();
				kabschFood[1] = c.getFragments().get(j).getAllResidues();
				t = Kabsch.calculateTransformation(kabschFood);
				temp += t.getRmsd();
			}
		}
		double result = temp / ((size - 1) * size / 2);
		return result;
	}

	// /**
	// * shortcut to retrieve a two-dimensional position from the
	// one-dimensional
	// * array.
	// *
	// * @param i
	// * the x-coordinate
	// * @param j
	// * the y-coordinate
	// * @return
	// */
	// private int getDistAt(int i, int j, int[] distance) {
	// if (j > i)
	// return distance[((j * (j - 1)) / 2) + i];
	// else
	// return distance[((i * (i - 1)) / 2) + j];
	// }

	public double calcInterClusterRMSD() {
		try {
			BufferedWriter w = new BufferedWriter(new FileWriter("/home/galicae/Documents/interCluster2"));
			BufferedWriter w1 = new BufferedWriter(new FileWriter("/home/galicae/Documents/interClusterGuide2"));
			double temp = 0;
			for (int i = 0; i < clusters.size(); i++) {
				w1.write(i + " " + clusters.get(i).getID() + "\n");
				for (int j = 0; j < clusters.size(); j++) {
					double cur = calcDistance(clusters.get(i), clusters.get(j));
					if(Double.isNaN(cur))
						cur=0;
//					if(cur < 0.4 && cur > 0.01)
//						System.err.println("bad cluster");
					temp += cur;
					w.write(cur + "\t");
				}
				w.write("\n");
			}
			double result = temp
					/ ((clusters.size() - 1) * clusters.size() / 2);
			w.close();
			w1.close();
			return result / 2;
		} catch (Exception e) {
			e.printStackTrace();
		}
		return -1;
	}

	private double calcDistance(FragmentCluster p, FragmentCluster q) {
		double[][][] kabschFood = new double[2][p.getFragmentLength()][3];
		kabschFood[0] = p.getCentroid().getAllResidues();
		kabschFood[1] = q.getCentroid().getAllResidues();
		Transformation t = Kabsch.calculateTransformation(kabschFood);
		return t.getRmsd();
	}
}
