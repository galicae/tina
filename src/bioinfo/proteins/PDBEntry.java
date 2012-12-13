/******************************************************************************
 * bioinfo.proteins.PDBEntry                                                  *
 *                                                                            *
 * This file is best read at line width 80 and tab width 4.                   *
 ******************************************************************************/
package bioinfo.proteins;

import java.io.Serializable;
import java.util.List;

import bioinfo.alignment.Alignable;

/**
 * @author gobi_4
 * @date November 24, 2012
 */
public class PDBEntry implements Alignable,Serializable{
	
	private final String id;
	private AminoAcid[] aminoAcids;
	private int length;
	
	/**
	 * Constructor for an the PDBEntry
	 * @param arg1 the id of the PDBEntry
	 */
	public PDBEntry(String arg1) {
		this.id = arg1;
	}
	
	/**
	 * Constructor for an the PDBEntry
	 * @param arg1 the id of the PDBEntry
	 * @param arg2 Array of AminoAcids
	 */
	public PDBEntry(String arg1, AminoAcid[] arg2) {
		this.id = arg1;
		this.aminoAcids = arg2;
		this.length = arg2.length;
	}
	
	/**
	 * Constructor for an the PDBEntry
	 * @param arg1 the id of the PDBEntry
	 * @param arg2 Array of AminoAcids
	 */
	public PDBEntry(String arg1, List<AminoAcid> arg2) {
		this.id = arg1;
		this.aminoAcids = arg2.toArray(new AminoAcid[0]);
		this.length = arg2.size();
	}
	
	/**
	 * @return the ID of the PDBEntry
	 */
	@Override
	public String getID() {
		return id;
	}
	
	/**
	 * Adds an AminoAcid to the entry. Usage discouraged because very slow.
	 * Better give the Constructor the full AminoAcid List.
	 * @deprecated
	 * @param arg
	 */
	public void addAminoAcid(AminoAcid arg) {
		AminoAcid[] result = new AminoAcid[length];
		for (int i = 0; i < length; i++) {
			result[i] = aminoAcids[i];
		}
		result[length] = arg;
		aminoAcids = result;
		length++;
	}
	
	public AminoAcid getAminoAcid(int n) {
		return aminoAcids[n];
	}

	@Override
	public int length() {
		return aminoAcids.length;
	}

	@Override
	public Object getComp(int n) {
		return getAminoAcid(n);
	}
	
	/**
	 * @return all AminoAcids as a (human readable) String
	 */
	public String getAtomSectionAsString(){
		int lineCounter = 0;
		String out = "";
		for(int i = 0; i != aminoAcids.length; i++){
			for(int j = 0; j != aminoAcids[i].getAtomNumber(); j++){
				out += aminoAcids[i].getAtom(j).toString(lineCounter++,i,aminoAcids[i].getName().getThreeLetterCode(),id.charAt(4)).trim()+"\n";
			}
		}
		return out;
	}
	
}
