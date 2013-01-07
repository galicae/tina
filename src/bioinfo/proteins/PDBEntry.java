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
public class PDBEntry implements Alignable, Serializable {

	private final String id;
	private final char chainID;
	private final int chainIDNum;
	private AminoAcid[] aminoAcids;
	private int length;

	/**
	 * Constructor for an the PDBEntry
	 * 
	 * @param arg1
	 *            the id of the PDBEntry
	 */
	public PDBEntry(String arg1) {
		if (arg1.length() == 4) {
			this.id = arg1;
			this.chainID = 'A';
			this.chainIDNum = 0;
		} else {
			this.id = arg1.substring(0, 4);
			this.chainID = arg1.charAt(4);
			this.chainIDNum = Integer.valueOf(arg1.substring(5, 7));
		}
	}

	/**
	 * Constructor for an the PDBEntry
	 * 
	 * @param arg1
	 *            the id of the PDBEntry
	 * @param arg2
	 *            Array of AminoAcids
	 */
	public PDBEntry(String arg1, AminoAcid[] arg2) {
		if (arg1.length() == 4) {
			this.id = arg1;
			this.chainID = 'A';
			this.chainIDNum = 0;
		} else {
			this.id = arg1.substring(0, 4);
			this.chainID = arg1.charAt(4);
			this.chainIDNum = Integer.valueOf(arg1.substring(5, 7));
		}
		this.aminoAcids = arg2;
		this.length = arg2.length;
	}

	/**
	 * Constructor for an the PDBEntry
	 * 
	 * @param arg1
	 *            the id of the PDBEntry
	 * @param arg2
	 *            Array of AminoAcids
	 */
	public PDBEntry(String arg1, List<AminoAcid> arg2) {
		if (arg1.length() == 4) {
			this.id = arg1;
			this.chainID = 'A';
			this.chainIDNum = 0;
		} else {
			this.id = arg1.substring(0, 4);
			this.chainID = arg1.charAt(4);
			this.chainIDNum = Integer.valueOf(arg1.substring(5));
			this.aminoAcids = arg2.toArray(new AminoAcid[0]);
			this.length = arg2.size();
		}
	}
	
	/**
	 * Constructor for an the PDBEntry
	 * 
	 * @param arg1
	 *            the id of the PDBEntry
	 * @param arg2
	 *            Array of AminoAcids
	 */
	public PDBEntry(String arg1, List<AminoAcid> arg2, boolean real) {
		if (arg1.length() == 4) {
			this.id = arg1;
			this.chainID = 'A';
			this.chainIDNum = 0;
		} else {
			this.id = arg1.substring(0, 4);
			this.chainID = arg1.charAt(4);
			this.chainIDNum = Integer.valueOf(arg1.substring(5));
			this.aminoAcids = arg2.toArray(new AminoAcid[0]);
			this.length = arg2.size();
		}
	}

	/**
	 * @return the ID of the PDBEntry
	 */
	@Override
	public String getID() {
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

	/**
	 * Adds an AminoAcid to the entry. Usage discouraged because very slow.
	 * Better give the Constructor the full AminoAcid List.
	 * 
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
	public String getAtomSectionAsString() {
		int lineCounter = 0;
		String out = "";
		for (int i = 0; i != aminoAcids.length; i++) {
			for (int j = 0; j != aminoAcids[i].getAtomNumber(); j++) {
				// TODO Some Entries have no 7-letter ID (1TIMA00) but four
				// letters (1TIM)
				// ==> make some Method for IDs without chainID!
				out += aminoAcids[i]
						.getAtom(j)
						.toString(lineCounter++, i,
								aminoAcids[i].getName().getThreeLetterCode(),
								chainID).trim()
						+ "\n";
			}
		}
		return out;
	}

	/**
	 * 
	 * @return all backbone atoms as Atom[]
	 * not existing atoms are written as null and have to be checked in all following methods
	 */
	public Atom[] getBackboneAtoms(){
		AminoAcid tmp;
		Atom[] backbone = new Atom[aminoAcids.length*4];
		Atom[] temp;
		for(int i = 0; i != aminoAcids.length; i++){
			tmp = aminoAcids[i];
			temp = tmp.getBackboneAtoms();
			for(int j = 0; j != 4; j++){
				backbone[i*4+j] = temp[j];
			}
		}
		return backbone;
	}
	
}
