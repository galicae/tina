package bioinfo.proteins;

import java.util.List;
import java.util.Locale;

import bioinfo.StructuralData;

/**
 * 
 * @author gobi_12_4
 * 
 */
public class DSSPEntry implements StructuralData{

	private final String id;
	private final char chainID;
	private final int chainIDNum;
	private int length;

	/**
	 * Sequence of AminoAcids
	 */
	private AminoAcidName[] names;
	/**
	 * Sequence of Secondary Structures
	 */
	private SecStructEight[] secondaryStructure;
	/**
	 * Sequence of solvent accessibility
	 */
	private int[] accesability;
	/**
	 * sequence of phi angles
	 */
	private double[] phi;
	/**
	 * sequence of psi angles
	 */
	private double[] psi;
	private double[][] caTrace;
	private int[] resIndex;

//	 public DSSPEntry(String arg1) {
//	 	if (arg1.length() == 4) {
//	 		this.id = arg1;
//	 		this.chainID = 'A';
//	 		this.chainIDNum = 0;
//	 	} else {
//	 		this.id = arg1.substring(0, 4);
//	 		this.chainID = arg1.charAt(4);
//	 		this.chainIDNum = Integer.valueOf(arg1.substring(5, 7));
//	 	}
//	 }

	/**
	 * Constructs a new DSSPEntry
	 * 
	 * @param id
	 *            the ID of the Protein
	 * @param names
	 *            the sequence of AminoAcid (-Names)
	 * @param resIndex
	 * @param secondaryStructure
	 * @param accessability
	 * @param phi
	 * @param psi
	 * @param caTrace
	 */
	public DSSPEntry(String id, AminoAcidName[] names, int[] resIndex,
			SecStructEight[] secondaryStructure, int[] accessability,
			double[] phi, double[] psi, double[][] caTrace) {
		if (id.length() == 4) {
			this.id = id;
			this.chainID = 'A';
			this.chainIDNum = 0;
		} else {
			this.id = id.substring(0, 4);
			this.chainID = id.charAt(4);
			this.chainIDNum = Integer.valueOf(id.substring(5, 7));
		}
		this.names = names;
		this.resIndex = resIndex;
		this.secondaryStructure = secondaryStructure;
		this.accesability = accessability;
		this.phi = phi;
		this.psi = psi;
		this.caTrace = caTrace;
		this.length = names.length;
	}

	/**
	 * convenience constructor
	 * 
	 * @param id
	 * @param names
	 * @param resIndex
	 * @param secondaryStructure
	 * @param accessability
	 * @param phi
	 * @param psi
	 * @param caTrace
	 */
	public DSSPEntry(String id, List<AminoAcidName> names,
			List<Integer> resIndex, List<SecStructEight> secondaryStructure,
			List<Integer> accessability, List<Double> phi, List<Double> psi,
			List<double[]> caTrace) {
		if (id.length() == 4) {
			this.id = id;
			this.chainID = 'A';
			this.chainIDNum = 0;
		} else {
			this.id = id.substring(0, 4);
			this.chainID = id.charAt(4);
			this.chainIDNum = Integer.valueOf(id.substring(5, 7));
		}
		this.length = names.size();
		this.names = names.toArray(new AminoAcidName[names.size()]);
		this.secondaryStructure = secondaryStructure
				.toArray(new SecStructEight[secondaryStructure.size()]);
		this.accesability = new int[accessability.size()];
		this.resIndex = new int[resIndex.size()];
		this.phi = new double[phi.size()];
		this.psi = new double[psi.size()];
		this.caTrace = new double[caTrace.size()][3];
		for (int i = 0; i != names.size(); i++) {
			this.resIndex[i] = resIndex.get(i);
			this.accesability[i] = accessability.get(i);
			this.phi[i] = phi.get(i);
			this.psi[i] = psi.get(i);
			this.caTrace[i] = caTrace.get(i);
		}
	}

	public int length() {
		return length;
	}

	public AminoAcidName[] getNames() {
		return names;
	}

	public int[] getResIndex() {
		return resIndex;
	}

	public SecStructEight[] getSecondaryStructure() {
		return secondaryStructure;
	}

	public int[] getAccesability() {
		return accesability;
	}

	public double[] getPhi() {
		return phi;
	}

	public double[] getPsi() {
		return psi;
	}

	public double[][] getCaTrace() {
		return caTrace;
	}

	/**
	 * return complete id
	 */
	public String getID(){
		return id+chainID+chainIDNum;
	}
	
	/**
	 * @return the ID
	 */
	public String getId() {
		return id;
	}

	/**
	 * @return the chainID
	 */
	public char getChainID() {
		return chainID;
	}

	/**
	 * @return the chainIDNum
	 */
	public int getChainIDNum() {
		return chainIDNum;
	}
	
	public AminoAcidName[] getSequenceData(){
		return this.getNames();
	}

	/**
	 * 
	 * @return String containing String representation of contents <br />
	 *         NOT DSSP File Format!
	 */
	public String asTable() {
		String out = "#AA\tSS\tACC\tPHI\tPSI\tCAX\tCAY\tCAZ\n";
		for (int i = 0; i != names.length; i++) {
			out += names[i].getOneLetterCode() + "\t";
			out += secondaryStructure[i].getCharRespres() + "\t";
			out += accesability[i] + "\t";
			out += String.format(Locale.US, "%.1f", phi[i]) + "\t";
			out += String.format(Locale.US, "%.1f", psi[i]) + "\t";
			out += String.format(Locale.US, "%.1f", caTrace[i][0]) + "\t";
			out += String.format(Locale.US, "%.1f", caTrace[i][1]) + "\t";
			out += String.format(Locale.US, "%.1f", caTrace[i][2]) + "\n";
		}
		return out;
	}

}
