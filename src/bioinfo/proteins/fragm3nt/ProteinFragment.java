package bioinfo.proteins.fragm3nt;

public class ProteinFragment {
	private String id;
	private String sequence;
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

	public ProteinFragment(String id, String seq, double[][] coordinates, int startIndex,
			int fragLength) {
		this.fragLength = fragLength;
		this.id = id;
		this.sequence = seq;
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
	
	public void setCoordinates(double[][] nCoord) {
		this.coordinates = nCoord;
	}
	public String toString() {
		String result = "";
		result += "Fragment " + id + ", residue length " + fragLength + "\n";
		for(int i = 0; i < coordinates.length; i++) {
			result += coordinates[i][0] + " " + coordinates[i][1] + " " + coordinates[i][2] + "\n";
		}
		return result;
	}
	
	public boolean equals(ProteinFragment other) {
		for (int i = 0; i < other.getFragmentLength(); i++) {
			for (int j = 0; j < 3; j++) {
				if (!isInEpsilon(coordinates[i][j], other.getAllResidues()[i][j]))
					return false;
			}
		}		
		return true;
	}
	
	private static final double epsilon = 0.0001d;

	private static boolean isInEpsilon(double a, double b) {
		return (a > (b - epsilon)) && (a < (b + epsilon));
	}

	public String getSequence() {
		return this.sequence;
	}
}
