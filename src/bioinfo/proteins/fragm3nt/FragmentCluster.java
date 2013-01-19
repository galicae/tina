package bioinfo.proteins.fragm3nt;

import java.util.LinkedList;

public class FragmentCluster {
	private LinkedList<ProteinFragment> fragments = new LinkedList<ProteinFragment>();
	private String name;
	private ProteinFragment centroid;
	private double[][] pssm;

	private static final double epsilon = 0.0001d;

	/**
	 * this function calculates the centroid of a cluster.
	 */
	public void calculateCentroid() {
		int clusterSize = fragments.size();
		double[][] newCentroid = new double[centroid.fragLength][3];
		ProteinFragment curFragment;
		double[] curResidue = new double[3];
		for (int j = 0; j < clusterSize; j++) {
			curFragment = fragments.get(j);
			for (int i = 0; i < curFragment.fragLength; i++) {
				curResidue = curFragment.getResidue(i);
				newCentroid[i][0] += curResidue[0];
				newCentroid[i][1] += curResidue[1];
				newCentroid[i][2] += curResidue[2];
			}
		}
		for (int i = 0; i < newCentroid.length; i++) {
			newCentroid[i][0] /= clusterSize * 1.;
			newCentroid[i][1] /= clusterSize * 1.;
			newCentroid[i][2] /= clusterSize * 1.;
		}

		centroid.setCoordinates(newCentroid);
	}

	/**
	 * Compares another cluster to this one by comparing their centroid
	 * coordinates. NOT an override of the equals(Object o) function.
	 * 
	 * @param f
	 *            the fragment cluster to compare with
	 * @return true if the coordinates are within an epsilon of each other;
	 *         false if not
	 */
	public boolean equals(FragmentCluster f) {
		if (f.getCentroid().getFragmentLength() != centroid.getFragmentLength())
			return false;
		for (int i = 0; i < centroid.getFragmentLength(); i++) {
			for (int j = 0; j < 3; j++) {
				if (!isInEpsilon(f.getCentroid().getAllResidues()[i][j],
						centroid.getAllResidues()[i][j]))
					return false;
			}
		}
		return true;
	}

	/**
	 * this function sets the cluster index of every fragment to -1 (not
	 * clustered) and removes them from the fragment list
	 */
	public void flush() {
		while (!fragments.isEmpty()) {
			fragments.remove();
		}
	}

	/**
	 * returns a string representation of the fragment cluster in pdb format.
	 * Each fragment in the cluster is wrapped as a pdb model.
	 */
	@Override
	public String toString() {
		int i = 0;
		StringBuilder result = new StringBuilder();
		try {
			for (ProteinFragment f : fragments) {
				result.append("REMARK 500 " + f.getSequence() + "\n");
				result.append("MODEL        " + ++i + "\n");
				result.append(f.toString());
				result.append("ENDMDL" + "\n");
			}
			return result.toString();
		} catch (Exception e) {
			e.printStackTrace();
		}
		return "problem printing";
		// the try-catch block owes its existence to InvocationException calls
		// from phantom functions. We still don't know what was wrong there, but
		// apparently we fixed it.
	}

	public int getSize() {
		return fragments.size();
	}

	public LinkedList<ProteinFragment> getFragments() {
		return fragments;
	}

	public double[][] getPssm() {
		return pssm;
	}

	public void setPssm(double[][] matrix) {
		pssm = matrix;
	}

	public String getID() {
		return name;
	}

	public ProteinFragment getCentroid() {
		return centroid;
	}

	public void setCentroid(ProteinFragment centroid) {
		this.centroid = centroid;
	}

	private static boolean isInEpsilon(double a, double b) {
		return (a > (b - epsilon)) && (a < (b + epsilon));
	}

	public void add(ProteinFragment f) {
		fragments.add(f);
	}
}