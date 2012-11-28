/******************************************************************************
 * bioinfo.proteins.AminoAcid                                                 *
 *                                                                            *
 * This file is best read at line width 80 and tab width 4.                   *
 ******************************************************************************/
package bioinfo.proteins;

import bioinfo.proteins.Atom.AtomType;

/**
 * @author gobi_4
 * @date November 24, 2012
 */
public class AminoAcid {

	/**
	 * enum with all the AminoAcidNames. Also contains methods to convert from
	 * one- to three-letter codes.
	 * @author gobi_4
	 */
	
	
	private final AminoAcidName name;
	private int numberOfAtoms;
	private Atom[] atoms;
	
	/**
	 * Constructor for the AminoAcid Class. Takes the AminoAcidName as arg.
	 * @param arg
	 */
	public AminoAcid(AminoAcidName arg) {
		this.name = arg;
	}
	
	/**
	 * Constructor for the AminoAcid Class. Takes the OneLetterCode of the
	 * AminoAcid and the Atom Array. 
	 * @param oneLetterCode
	 * @param atoms
	 */
	public AminoAcid(AminoAcidName arg, Atom[] atoms) {
		this.name = arg;
		this.atoms = atoms;
		this.numberOfAtoms = atoms.length;
	}
	
	/**
	 * Constructor for the AminoAcid Class. Takes the OneLetterCode of the
	 * AminoAcid as arg.
	 * @param oneLetterCode
	 */
	public AminoAcid(String oneLetterCode) {
		this.name = AminoAcidName.getAAFromOLC(oneLetterCode);
	}
	
	/**
	 * Constructor for the AminoAcid Class. Takes the OneLetterCode of the
	 * AminoAcid and the Atom Array. 
	 * @param oneLetterCode
	 * @param atoms
	 */
	public AminoAcid(String oneLetterCode, Atom[] atoms) {
		this.name = AminoAcidName.getAAFromOLC(oneLetterCode);
		this.atoms = atoms;
		this.numberOfAtoms = atoms.length;
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
	 * @param n
	 * @return the n-th Atom
	 */
	public Atom getAtom(int n) {
		if (n >= numberOfAtoms) {
			System.err.println("This Atom doesn't exist: "+n);
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
		System.err.println("The given Atom doesn't exist: " + atomtype);
		return null;
	}
	
	/**
	 * adds an atom to the aminoAcid. Usage discouraged as very slow!
	 * @deprecated
	 * @param atom
	 */
	public void addAtom(Atom atom) {
		Atom[] result = new Atom[numberOfAtoms+1];
		for (int i = 0; i < numberOfAtoms; i++) {
			result[i] = atoms[i];
		}
		result[numberOfAtoms] = atom;
		atoms = result;
		numberOfAtoms++;
	}
	
}
