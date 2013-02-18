package bioinfo.proteins;

/**
 * enum with all the AminoAcidNames. Also contains methods to convert from
 * one- to three-letter codes.
 * @author gobi_4
 */	
public enum AminoAcidName {
	A("A", "ALA", "Alanine", 0), 
	R("R", "ARG", "Arginine", 17), 
	N("N", "ASN","Asparagine", 13), 
	D("D", "ASP", "Aspartic acid", 3), 
	C("C", "CYS","Cysteine", 2), 
	E("E", "GLU", "Glutamic acid", 4), 
	Q("Q", "GLN","Glutamine", 16), 
	G("G", "GLY", "Gylcine", 6), 
	H("H", "HIS", "Histidine", 7), 
	I("I", "ILE", "Isoleucine", 8), 
	L("L", "LEU", "Leucine", 11), 
	K("K", "LYS","Lysine", 10),
	M("M", "MET", "Methionine", 12), 
	F("F", "PHE","Phenylalanine", 5), 
	P("P", "PRO", "Proline", 15), 
	S("S", "SER", "Serine", 18),
	T("T", "THR", "Threonine", 19), 
	W("W", "TRP", "Tryptophane", 22), 
	Y("Y", "TYR", "Tyrosine", 24), 
	V("V", "VAL", "Valine", 21),
	U("U", "UKN", "Unknown", 20);

	private final String oneLetterCode;
	private final String threeLetterCode;
	private final String longName;
	private final int number;

	/**
	 * Constructor
	 * 
	 * @param oneLetterCode
	 * @param threeLetterCode
	 * @param longName
	 */
	AminoAcidName(String oneLetterCode, String threeLetterCode, String longName, int number) {
		this.oneLetterCode = oneLetterCode;
		this.threeLetterCode = threeLetterCode;
		this.longName = longName;
		this.number = number;
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
	
	public int getNumber() {
		return number;
	}

	public static AminoAcidName getAAFromOLC(String oneLetterCode) {
		for (AminoAcidName aa : AminoAcidName.values()) {
			if (aa.getOneLetterCode().equals(oneLetterCode))
				return aa;
		}
		return U;
	}
	
	public static AminoAcidName getAAFromOLC(char oneLetterCode) {
		for (AminoAcidName aa : AminoAcidName.values()) {
			if (aa.getOneLetterCode().equals(String.valueOf(oneLetterCode)))
				return aa;
		}
		return U;
	}

	public static AminoAcidName getAAFromTLC(String threeLetterCode) {
		for (AminoAcidName aa : AminoAcidName.values()) {
			if (aa.getThreeLetterCode().equals(threeLetterCode))
				return aa;
		}
		return U;
	}
	
	public static AminoAcidName getAAFromNumber(int number) {
		for (AminoAcidName aa : AminoAcidName.values()) {
			if (aa.getNumber() == number)
				return aa;
		}
		return U;
	}
	
	public static int getNumberFromTLC(String threeLetterCode) {
		for (AminoAcidName aa : AminoAcidName.values()) {
			if (aa.getThreeLetterCode().equals(threeLetterCode))
				return aa.getNumber();
		}
		return U.getNumber();
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

	public static String getTLCFromNumber(int number) {
		for (AminoAcidName aa : AminoAcidName.values()) {
			if (aa.getNumber() == number)
				return aa.threeLetterCode;
		}
		return null;
	}
	
}
