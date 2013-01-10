package bioinfo.proteins.fragm3nt;

import bioinfo.proteins.Atom;
import bioinfo.proteins.AtomType;

public class ProteinFragment {
	private String id;
	private String sequence;
	private double[][] coordinates;
	private Atom[] atoms;
	private int startIndex;
	public final int fragLength;
	private boolean visited = false;
	private boolean noise = false;
	private int clustered = -1; // this should save the position of the cluster
								// to which the fragment belongs. It will save a
								// lot of time during the initialization of the
								// clusters.

	public ProteinFragment(String id, Atom[] atoms, int startIndex,
			int fragLength) {
		this.atoms = atoms;
		this.fragLength = fragLength;
		this.id = id;
		this.coordinates = new double[atoms.length][3];
		this.startIndex = startIndex;
		for(int i = 0; i < atoms.length; i++) {
			coordinates[i] = atoms[i].getPosition();
		}
	}

	// centroid constructor
	public ProteinFragment(String id, double[][] atoms, int startIndex,
			int fragLength) {
		this.fragLength = fragLength;
		this.id = id;
		this.coordinates = new double[atoms.length][3];
		this.coordinates = atoms;
		this.startIndex = startIndex;
	}
	
	public ProteinFragment(String id, String seq, Atom[] atoms,
			int startIndex, int fragLength) {
		this.atoms = atoms;
		this.fragLength = fragLength;
		this.id = id;
		this.sequence = seq;
		this.coordinates = new double[atoms.length][3];
		this.startIndex = startIndex;
		for(int i = 0; i < atoms.length; i++) {
			coordinates[i] = atoms[i].getPosition();
		}
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
		for(int i = 0; i < nCoord.length; i++) {
			atoms[i].setPosition(nCoord[i]);
		}
	}

	public String toString() {
		String result = "";
		result += "Fragment " + id + ", residue length " + fragLength + "\n";
		for (int i = 0; i < coordinates.length; i++) {
			result += coordinates[i][0] + " " + coordinates[i][1] + " "
					+ coordinates[i][2] + "\n";
		}
		return result;
	}
	
	public String getCanonicalRepresentation() {
		StringBuilder sb = new StringBuilder();
		Atom a = new Atom(AtomType.C, new double[2]);
		for(int i = 0; i < atoms.length; i++) {
			a = atoms[i];
			sb.append(a.toString(i, (i/4), String.valueOf(sequence.charAt(i/4)), 'A') + "\n");
		}
		
//		int i = atoms.length;
//		String term = atoms[atoms.length - 1].toString(i, (i-1)/4, String.valueOf(sequence.charAt((i-1)/4)), 'A');
//		term = term.replace("ATOM", "TER ");
//		term = term.substring(0, 11) + "      " + term.substring(17, 25);
//		sb.append(term + "\n");
		return sb.toString();
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

	public void setNoise(boolean noise) {
		this.noise = noise;
	}
	
	public Atom[] getAtoms() {
		return this.atoms;
	}
	
	public void correctCoordinates(double[] correction) {
		for(int i = 0; i < coordinates.length; i++) {
			for(int j = 0; j < coordinates[0].length; j++) {
				coordinates[i][j] -= correction[j];
			}
		}
	}
}
