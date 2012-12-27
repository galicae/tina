package bioinfo.proteins.framg3nt;

import java.util.LinkedList;

public class FragmentCluster {
	private LinkedList<ProteinFragment> fragments = new LinkedList<ProteinFragment>();
	private String name;
	private ProteinFragment centroid;

	public void calculateCentroid() {
		int fragLength = centroid.getFragmentLength();
		double[][] newCentroid = new double[fragLength][3];
		ProteinFragment curFragment;
		double[] curResidue = new double[3];
		while (!fragments.isEmpty()) {
			curFragment = fragments.pop();
			for (int i = 0; i < fragLength; i++) {
				curResidue = curFragment.getResidue(i);
				newCentroid[i][0] += curResidue[0];
				newCentroid[i][1] += curResidue[1];
				newCentroid[i][2] += curResidue[2];
			}
		}
		for (int i = 0; i < fragLength; i++) {
			newCentroid[i][0] /= fragLength * 1.;
			newCentroid[i][1] /= fragLength * 1.;
			newCentroid[i][2] /= fragLength * 1.;
		}

		centroid = new ProteinFragment(centroid.getID(), newCentroid,
				centroid.getStartIndex(), fragLength);
	}

	public ProteinFragment getCentroid() {
		return centroid;
	}
	
	public void setCentroid(ProteinFragment centroid) {
		this.centroid = centroid;
	}

	@Override
	public boolean equals(Object o) {
		FragmentCluster f = (FragmentCluster) o;
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

	private static final double epsilon = 0.0001d;

	private static boolean isInEpsilon(double a, double b) {
		return (a > (b - epsilon)) && (a < (b + epsilon));
	}
}
