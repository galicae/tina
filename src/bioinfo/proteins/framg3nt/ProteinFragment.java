package bioinfo.proteins.framg3nt;

public class ProteinFragment {
	private String id;
	private double[][] coordinates;
	private int startIndex;
	public final int fragLength;
	private int clustered = -1; // this should save the position of the cluster
								// to which the fragment belongs. It will save a
								// lot of time during the initialization of the
								// clusters.

	public ProteinFragment(String id, double[][] coordinates, int startIndex,
			int fragLength) {
		this.fragLength = fragLength;
		this.id = id;
		this.coordinates = new double[fragLength][3];
		this.coordinates = coordinates;
		this.startIndex = startIndex;
	}

	public double[] getResidue(int i) {
		return coordinates[i];
	}

	public double[][] getAllResidues() {
		return coordinates;
	}

	public int getFragmentLength() {
		return fragLength;
	}

	public String getID() {
		return id;
	}

	public int getStartIndex() {
		return startIndex;
	}

	public void setClusterIndex(int i) {
		clustered = i;
	}

	public int getClusterIndex() {
		return clustered;
	}
}
