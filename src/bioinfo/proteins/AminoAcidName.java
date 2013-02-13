package bioinfo.proteins;

/**
 * enum with all the AminoAcidNames. Also contains methods to convert from
 * one- to three-letter codes.
 * @author gobi_4
 */	
public enum AminoAcidName {
	A("A", "ALA", "Alanine"), 
	R("R", "ARG", "Arginine"), 
	N("N", "ASN","Asparagine"), 
	D("D", "ASP", "Aspartic acid"), 
	C("C", "CYS","Cysteine"), 
	E("E", "GLU", "Glutamic acid"), 
	Q("Q", "GLN","Glutamine"), 
	G("G", "GLY", "Gylcine"), 
	H("H", "HIS", "Histidine"), 
	I("I", "ILE", "Isoleucine"), 
	L("L", "LEU", "Leucine"), 
	K("K", "LYS","Lysine"),
	M("M", "MET", "Methionine"), 
	F("F", "PHE","Phenylalanine"), 
	P("P", "PRO", "Proline"), 
	S("S", "SER", "Serine"),
	T("T", "THR", "Threonine"), 
	W("W", "TRP", "Tryptophane"), 
	Y("Y", "TYR", "Tyrosine"), 
	V("V", "VAL", "Valine"),
	U("U", "UKN", "Unknown");

	private final String oneLetterCode;
	private final String threeLetterCode;
	private final String longName;

	/**
	 * Constructor
	 * 
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
	
	public static AminoAcidName getAAFromOLC(char oneLetterCode) {
		for (AminoAcidName aa : AminoAcidName.values()) {
			if (aa.getOneLetterCode().equals(String.valueOf(oneLetterCode)))
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
