package bioinfo;

import bioinfo.proteins.AminoAcidName;

public interface StructuralData {
	
	/**
	 * @return length of "main" data in this object
	 * eg amino acids, residues, contacts etc
	 */
	int length();
	
	/**
	 * 
	 * @return String representing id of this object
	 * in normal / most cases it will be the PDB id
	 */
	String getID();
	
	/**
	 * returns sequence Data as AminoAcidName[]
	 */
	AminoAcidName[] getSequenceData();

}
