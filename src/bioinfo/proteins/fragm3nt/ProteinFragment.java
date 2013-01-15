package bioinfo.proteins.fragm3nt;

import bioinfo.proteins.Atom;
import bioinfo.proteins.AtomType;

public class ProteinFragment {
	private String id;
	private String sequence;
	private double[][] coordinates;
	private int startIndex;
	public final int fragLength;
	private boolean visited = false;
	private boolean noise = false;
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

	public ProteinFragment(String id, String seq, double[][] coordinates,
			int startIndex, int fragLength) {
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

	/**
	 * override the normal toString() so that we get PDB Atom lines with correct
	 * spacing and amino acid names
	 */
	public String toString() {
		StringBuilder result = new StringBuilder();
		for (int i = 0; i < coordinates.length; i++) {
			AtomType type = AtomType.C;
			switch (i % 4) {
			case 0:
				type = AtomType.N;
				break;
			case 1:
				type = AtomType.CA;
				break;
			case 2:
				type = AtomType.O;
				break;
			case 3:
				type = AtomType.C;
				break;
			}
			Atom cur = new Atom(type, coordinates[i]);
			result.append(cur.toString(i, (i / 4),
					String.valueOf(sequence.charAt(i / 4)), 'A')
					+ "\n");
		}
		return result.toString();
	}

	/**
	 * parameterized toString function, for when we only need part of the atoms
	 * 
	 * @param start
	 *            the point from which to start in the array
	 * @return PDB atom lines with the coordinates and the correct amino acid
	 *         names and spacing
	 */
	public String toString(int start, int position) {
		if (start < 0)
			start = 0;
		StringBuilder result = new StringBuilder();
		for (int i = start * 4; i < coordinates.length; i++) {
			int corr = 0;
			AtomType type = AtomType.C;
			switch (i % 4) {
			case 0:
				type = AtomType.N;
				break;
			case 1:
				type = AtomType.CA;
				break;
			case 2:
				type = AtomType.O;
				break;
			case 3:
				type = AtomType.C;
				break;
			}
			Atom cur = new Atom(type, coordinates[i]);
			result.append(cur.toString(corr + i, ( (corr + i) / 4),
					String.valueOf(sequence.charAt(i / 4)), 'A')
					+ "\n");
		}
		return result.toString();
	}

	public boolean equals(ProteinFragment other) {
		for (int i = 0; i < other.getFragmentLength(); i++) {
			for (int j = 0; j < 3; j++) {
				if (!isInEpsilon(coordinates[i][j],
						other.getAllResidues()[i][j]))
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

	public boolean isVisited() {
		return visited;
	}

	public void setVisited(boolean vis) {
		visited = vis;
	}

	public boolean isNoise() {
		return noise;
	}

	public void setSequence(String seq) {
		sequence = seq;
	}

	public void setStartIndex(int i) {
		startIndex = i;
	}

	public void setNoise(boolean noise) {
		this.noise = noise;
	}

	@Override
	public ProteinFragment clone() {
		ProteinFragment result = new ProteinFragment(id, sequence, coordinates,
				startIndex, fragLength);
		return result;
	}
	
	public void correctCoordinates() {
		String coord = String.valueOf(coordinates[0][0]);
		if(coord.length() > 7) {
			double corr = coordinates[0][0];
			for(int i = 0; i < coordinates.length; i++) {
				for(int j = 0; j < coordinates[0].length; j++) {
					coordinates[i][j] -= corr;
				}
			}
		}
		
	}
}
