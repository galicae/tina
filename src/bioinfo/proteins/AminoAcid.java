/******************************************************************************
 * bioinfo.proteins.AminoAcid                                                 *
 *                                                                            *
 * This file is best read at line width 80 and tab width 4.                   *
 ******************************************************************************/
package bioinfo.proteins;

import java.io.Serializable;
//import java.util.ArrayList;
//import java.util.List;

/**
 * @author gobi_4
 * @date November 24, 2012
 */
public class AminoAcid implements Serializable {

	/**
	 * 
	 */
	private static final long serialVersionUID = 6808776190281449279L;

	/**
	 * name of the AminoAcid
	 */
	private final AminoAcidName name;
	/**
	 * Residue Index: index of the AminoAcid in the protein (for PDB output)
	 */
	private int resIndex;
	/**
	 * number of atoms currently in the AminoAcid
	 */
	private int numberOfAtoms;
	/**
	 * Array of Atoms currently assigned to this AminoAcid
	 */
	private Atom[] atoms;

	/**
	 * Standard constructor for an AminoAcid
	 * 
	 * @param aminoAcidName
	 * @param oneLetterCode
	 * @param atoms
	 */
	public AminoAcid(AminoAcidName aminoAcidName, int resIndex, Atom[] atoms) {
		this.name = aminoAcidName;
		this.resIndex = resIndex;
		setAtoms(atoms);
	}

	/**
	 * Constructor for an AminoAcid
	 * 
	 * @param aminoAcidName
	 * @param resIndex
	 */
	public AminoAcid(AminoAcidName aminoAcidName, int resIndex) {
		this(aminoAcidName, resIndex, null);
	}

	/**
	 * Constructor for an AminoAcid
	 * 
	 * @param oneLetterCode
	 * @param resIndex
	 */
	public AminoAcid(String oneLetterCode, int resIndex) {
		this(AminoAcidName.getAAFromOLC(oneLetterCode), resIndex);
	}

	/**
	 * Constructor for an AminoAcid
	 * 
	 * @param oneLetterCode
	 * @param resIndex
	 */
	public AminoAcid(char oneLetterCode, int resIndex) {
		this(AminoAcidName.getAAFromOLC(oneLetterCode), resIndex);
	}

	/**
	 * Constructor for an AminoAcid
	 * 
	 * @param oneLetterCode
	 * @param resIndex
	 * @param atoms
	 */
	public AminoAcid(String oneLetterCode, int resIndex, Atom[] atoms) {
		this(AminoAcidName.getAAFromOLC(oneLetterCode), resIndex, atoms);
	}

	/**
	 * Constructor for an AminoAcid
	 * 
	 * @param olc
	 * @param resIndex
	 * @param atoms
	 */
	public AminoAcid(char olc, int resIndex, Atom[] atoms) {
		this(AminoAcidName.getAAFromOLC(olc), resIndex, atoms);
	}

	/**
	 * 
	 * @return the AminoAcidName
	 */
	public AminoAcidName getName() {
		return name;
	}

	/**
	 * 
	 * @return the AminoAcidResidueIndex
	 */
	public int getResIndex() {
		return resIndex;
	}

	/**
	 * sets the AtomArray and updates this.length
	 * 
	 * @param atoms
	 */
	public void setAtoms(Atom[] atoms) {
		this.atoms = atoms;
		if (atoms != null) {
			this.numberOfAtoms = this.atoms.length;
		} else
			this.numberOfAtoms = 0;
	}

	/**
	 * 
	 * @param n
	 * @return the n-th Atom
	 */
	public Atom getAtom(int n) {
		if (n >= numberOfAtoms) {
			System.err.println("This Atom doesn't exist: " + n);
			return null;
		}
		return atoms[n];
	}

	/**
	 * 
	 * @param atomtype
	 * @return the atom of the given AtomType
	 */
	public Atom getAtomByType(AtomType atomtype) {
		for (int i = 0; i < numberOfAtoms; i++) {
			if (atoms[i].getType() == atomtype)
				return atoms[i];
		}
		// System.err.println("The given Atom doesn't exist: " + atomtype);
		return null;
	}

	/**
	 * adds an atom to the aminoAcid. Usage discouraged as very slow!
	 * 
	 * @deprecated
	 * @param atom
	 */
	public void addAtom(Atom atom) {
		Atom[] result = new Atom[numberOfAtoms + 1];
		for (int i = 0; i < numberOfAtoms; i++) {
			result[i] = atoms[i];
		}
		result[numberOfAtoms] = atom;
		atoms = result;
		numberOfAtoms++;
	}

	/**
	 * returns One-letter String representation of AminoAcid
	 */
	@Override
	public String toString() {
		return name.getOneLetterCode();
	}

	// 1 2 3 4 5 6 7 8
	// 12345678901234567890123456789012345678901234567890123456789012345678901234567890
	// MODEL 1
	// ATOM 1 N ALA A 1 11.104 6.134 -6.504 1.00 0.00 N
	// ATOM
	/**
	 * should return all atoms of an amino acid in PDB format
	 */
	public String toPDBLineString(int startIndex, char chain) {
		StringBuilder result = new StringBuilder();
		for (int i = 0; i < numberOfAtoms; i++) {
			result.append(atoms[i].toString((startIndex + i), resIndex,
					name.getOneLetterCode(), chain));
		}
		return result.toString();
	}

	/**
	 * 
	 * @return number of Atoms in this AminoAcid
	 */
	public int getAtomNumber() {
		return numberOfAtoms;
	}

	public AminoAcid copy() {
		return new AminoAcid(name, resIndex, atoms.clone());
	}

	/**
	 * 
	 * @return Atom[] containing all existing backbone atoms (best case all 4)
	 *         non-existing atoms are replaced with null and have to be checked
	 *         in all following methods
	 */
	public Atom[] getBackboneAtoms() {
		Atom[] result = new Atom[4];
		for (int i = 0; i != atoms.length; i++) {
			if (atoms[i] != null) {
				switch (atoms[i].getType()) {
				case C:
					result[2] = atoms[i];
					break;
				case CA:
					result[1] = atoms[i];
					break;
				case N:
					result[0] = atoms[i];
					break;
				case O:
					result[3] = atoms[i];
					break;
				default:
					break;
				}
			}
		}
		return result;
	}

	/**
	 * 
	 * @return
	 */
	public Atom[] getAtoms() {
		return this.atoms;
	}
}
