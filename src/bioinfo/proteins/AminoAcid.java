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
	public enum AminoAcidName {
		A ("A", "ALA", "Alanine"),
		R ("R", "ARG", "Arginine"),
		N ("N", "ASN", "Asparagine"),
		D ("D", "ASP", "Aspartic acid"),
		C ("C", "CYS", "Cysteine"),
		E ("E", "GLU", "Glutamic acid"),
		Q ("Q", "GLN", "Glutamine"),
		G ("G", "GLY", "Gylcine"),
		H ("H", "HIS", "Histidine"),
		I ("I", "ILE", "Isoleucine"),
		L ("L", "LEU", "Leucine"),
		K ("K", "LYS", "Lysine"),
		M ("M", "MET", "Methionine"),
		F ("F","PHE","Phenylalanine"),
		P ("P","PRO","Proline"),
		S ("S","SER","Serine"),
		T ("T","THR","Threonine"),
		W ("W","TRP","Tryptophane"),
		Y ("Y","TYR","Tyrosine"),
		V ("V","VAL","Valine");
		
		private final String oneLetterCode;
		private final String threeLetterCode;
		private final String longName;
		/**
		 * Constructor
		 * @param oneLetterCode
		 * @param threeLetterCode
		 * @param longName
		 */
		AminoAcidName(String oneLetterCode, String threeLetterCode, String longName) {
			this.oneLetterCode = oneLetterCode;
			this.threeLetterCode = threeLetterCode;
			this.longName = longName;
		}
		
		public String getOneLetterCode() {
			return oneLetterCode;
		}
		public String getThreeLetterCode() {
			return threeLetterCode;
		}
		public String getLongName() {
			return longName;
		}
		
		public static AminoAcidName getAAFromOLC(String oneLetterCode) {
			for (AminoAcidName aa : AminoAcidName.values()) {
				if (aa.getOneLetterCode().equals(oneLetterCode))
					return aa;
			}
			return null;
		}
		public static AminoAcidName getAAFromTLC(String threeLetterCode) {
			for (AminoAcidName aa : AminoAcidName.values()) {
				if (aa.getThreeLetterCode().equals(threeLetterCode))
					return aa;
			}
			return null;
		}
		
		public static String getThreeLetterCode(String oneLetterCode) {
			for (AminoAcidName aa : AminoAcidName.values()) {
				if (aa.getOneLetterCode().equals(oneLetterCode))
					return aa.threeLetterCode;
			}
			return null;
		}
		public static String getOneLetterCode(String threeLetterCode) {
			for (AminoAcidName aa : AminoAcidName.values()) {
				if (aa.getThreeLetterCode().equals(threeLetterCode))
					return aa.oneLetterCode;
			}
			return null;
		}
	}
	
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
	 * adds an atom to the aminoAcid. Usage discourages as very slow!
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
