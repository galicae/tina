package bioinfo.proteins.fragm3nt;

import bioinfo.proteins.AminoAcid;
import bioinfo.proteins.Atom;
import bioinfo.proteins.AtomType;
import bioinfo.proteins.PDBEntry;

public class ProteinFragment {
	private String id;
	private String sequence;
	private double[][] coordinates;
	private Atom[] atoms;
	public final int fragLength;
	private boolean visited = false;
	private boolean noise = false;
	private int clustered = -1; // this should save the position of the cluster
								// to which the fragment belongs. It will save a
								// lot of time during the initialization of the
								// clusters.

	public ProteinFragment(String id, Atom[] atoms, int fragLength) {
		this.fragLength = fragLength;
		this.id = id;
		this.atoms = atoms;
		coordinates = new double[atoms.length][3];
		for (int i = 0; i < atoms.length; i++) {
			coordinates[i] = atoms[i].getPosition();
		}
	}

	public ProteinFragment(String id, String seq, Atom[] atoms, int fragLength) {
		this.fragLength = fragLength;
		this.id = id;
		this.sequence = seq;
		this.atoms = atoms;
		coordinates = new double[atoms.length][3];
		for (int i = 0; i < atoms.length; i++) {
			coordinates[i] = atoms[i].getPosition();
		}
	}
	
	public ProteinFragment(String id, String seq, double[][] newCentroid, int fragLength) {
		this.fragLength = fragLength;
		this.id = id;
		this.sequence = seq;
		coordinates = new double[newCentroid.length][3];
		atoms = new Atom[newCentroid.length];
		for (int i = 0; i < newCentroid.length; i++) {
			atoms[i] = new Atom(AtomType.CA, newCentroid[i]);
			coordinates[i] = newCentroid[i];
		}
	}

	public ProteinFragment(String id, double[][] newCentroid, int fragLength) {
		this.fragLength = fragLength;
		this.id = id;
		coordinates = new double[newCentroid.length][3];
		atoms = new Atom[newCentroid.length];
		for (int i = 0; i < newCentroid.length; i++) {
			atoms[i] = new Atom(AtomType.CA, newCentroid[i]);
			coordinates[i] = newCentroid[i];
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

	public void setClusterIndex(int i) {
		clustered = i;
	}

	public int getClusterIndex() {
		return clustered;
	}

	public void setCoordinates(double[][] nCoord) {
		// this.coordinates = nCoord;
		for (int i = 0; i < nCoord.length; i++) {
			coordinates[i] = nCoord[i];
			atoms[i].setPosition(nCoord[i]);
		}
	}

	/**
	 * override the normal toString() so that we get PDB Atom lines with correct
	 * spacing and amino acid names
	 */
	public String toString() {
		StringBuilder result = new StringBuilder();
		for (int i = 0; i < coordinates.length; i++) {
			result.append(atoms[i].toString(i, i,
					String.valueOf(sequence.charAt(i)), 'A')
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
		for (int i = start; i < coordinates.length; i++) {
			int ind = position + i - start;
			result.append(atoms[i].toString(ind, ind,
					String.valueOf(sequence.charAt(i)), 'A')
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

	public void setNoise(boolean noise) {
		this.noise = noise;
	}

	@Override
	public ProteinFragment clone() {
		String id = this.id;
		String sequence = this.sequence;
		Atom[] atoms = new Atom[this.atoms.length];
		for (int i = 0; i < atoms.length; i++) {
			atoms[i] = new Atom(this.atoms[i].getType(), coordinates[i]);
		}
		ProteinFragment result = new ProteinFragment(id, sequence, atoms,
				fragLength);
		result.setClusterIndex(clustered);
		return result;
	}

	public void translateCoordinates(double[] correct) {
		for (int i = 0; i < coordinates.length; i++) {
			for (int j = 0; j < coordinates[0].length; j++) {
				coordinates[i][j] += correct[j];
				atoms[i].setPosition(coordinates[i]);
			}
		}
	}

	public void append(double[][] coord, String seq) {
		double[][] newCoordinates = new double[coord.length
				+ coordinates.length][3];
		Atom[] newAtoms = new Atom[coord.length + coordinates.length];
		int i = 0;

		for (i = 0; i < coordinates.length; i++) {
			newCoordinates[i] = coordinates[i];
			newAtoms[i] = atoms[i];
		}
		for (int j = i; j < newCoordinates.length; j++) {
			newCoordinates[j] = coord[j - i];
			newAtoms[j] = new Atom(AtomType.CA, coord[j - i]);
		}

		this.coordinates = new double[coord.length + coordinates.length][3];
		this.atoms = new Atom[coord.length + coordinates.length];
		for (int j = 0; j < coordinates.length; j++) {
			atoms[j] = newAtoms[j];
			for (int k = 0; k < 3; k++) {
				coordinates[j][k] = newCoordinates[j][k];
			}
		}
		this.sequence += seq;
	}

	public Atom[] getAtoms() {
		return this.atoms;
	}
}
