package bioinfo.proteins.fragm3nt;

import java.text.DecimalFormat;
import java.util.LinkedList;
import java.util.Locale;

public class FragmentCluster {
	private LinkedList<ProteinFragment> fragments = new LinkedList<ProteinFragment>();
	private String name;
	private ProteinFragment centroid;
	public double[][] pssm;

	public void calculateCentroid() {
		int actLength = centroid.getAllResidues().length;
		double[][] newCentroid = new double[actLength][3];
		ProteinFragment curFragment;
		double[] curResidue = new double[3];
		while (!fragments.isEmpty()) {
			curFragment = fragments.pop();
			for (int i = 0; i < actLength; i++) {
				curResidue = curFragment.getResidue(i);
				newCentroid[i][0] += curResidue[0];
				newCentroid[i][1] += curResidue[1];
				newCentroid[i][2] += curResidue[2];
			}
		}
		for (int i = 0; i < actLength; i++) {
			newCentroid[i][0] /= actLength * 1.;
			newCentroid[i][1] /= actLength * 1.;
			newCentroid[i][2] /= actLength * 1.;
		}

		centroid = new ProteinFragment(centroid.getID(), newCentroid,
				centroid.getStartIndex(), actLength);
	}

	public ProteinFragment getCentroid() {
		return centroid;
	}

	public void setCentroid(ProteinFragment centroid) {
		this.centroid = centroid;
	}

	public boolean equals(FragmentCluster o) {
		FragmentCluster f = o;
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

	public void add(ProteinFragment f) {
		fragments.add(f);
	}

	public void flush() {
		while (!fragments.isEmpty()) {
			// reset everything
			fragments.getFirst().setClusterIndex(-1);
			// then remove it
			fragments.remove();
		}
	}

	private static final double epsilon = 0.0001d;

	private static boolean isInEpsilon(double a, double b) {
		return (a > (b - epsilon)) && (a < (b + epsilon));
	}

	@Override
	public String toString() {
		StringBuilder result = new StringBuilder();
		ProteinFragment f = new ProteinFragment(null, new double[1][1], 0, 0);
		for(int i = 0; i < fragments.size(); i++) {
			f = fragments.get(i);
			result.append("REMARK " + f.getSequence() + "\n");
			result.append("MODEL        " + (i+1) + "\n");
			result.append(f.getCanonicalRepresentation());
			result.append("ENDMDL\n");
		}
		return result.toString();
	}

	public int getSize() {
		return fragments.size();
	}
	
	public LinkedList<ProteinFragment> getFragments() {
		return this.fragments;
	}
	
	public void setPssm(double[][] matrix) {
		pssm = matrix;
	}
}
