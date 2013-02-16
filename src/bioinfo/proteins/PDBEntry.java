/******************************************************************************
 * bioinfo.proteins.PDBEntry                                                  *
 *                                                                            *
 * Provides access to PDB Entries and also is our format for AminoAcid chain  *
 * structures                                                                 * 
 *                                                                            *
 * This file is best read at line width 80 and tab width 4.                   *
 ******************************************************************************/
package bioinfo.proteins;

import java.io.Serializable;
import java.util.List;

import bioinfo.alignment.Alignable;

/**
 * @author gobi_4
 * @lastchange 2013-02-16
 */
public class PDBEntry implements Alignable, Serializable {

	/**
	 * needed for class check. for serialization. please increase in case of
	 * editing this file
	 */
	private static final long serialVersionUID = 2L;
	
	private final String id;
	private final char chainID;
	private final int chainIDNum;
	private AminoAcid[] aminoAcids;
	private int length;

	/**
	 * Constructor for an a PDBEntry with an Array of AminoAcids
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
	 * Constructor for an the PDBEntry without any AminoAcids
	 * 
	 * @param arg1
	 *            the id of the PDBEntry
	 */
	public PDBEntry(String arg1) {
		this(arg1, new AminoAcid[0]);
	}
	
	/**
	 * Constructor for an the PDBEntry with a List of AminoAcids
	 * 
	 * @param arg1
	 *            the id of the PDBEntry
	 * @param arg2
	 *            Array of AminoAcids
	 */
	public PDBEntry(String arg1, List<AminoAcid> arg2) {
		this(arg1, arg2.toArray(new AminoAcid[0]));
	}

	/**
	 * Constructor for an an PDBEntry
	 * 
	 * @param arg1
	 *            the id of the PDBEntry
	 * @param arg2
	 *            Array of AminoAcids
	 * @param real
	 *            ?????
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

	/**
	 * returns the n-th AminoAcid
	 * @param n
	 * @return the n-th AminoAcid
	 */
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
	 * returns the AmioAcids as a (human readable) String like in a PDB file
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
				out += aminoAcids[i].getAtom(j).toString(lineCounter++, i, aminoAcids[i].getName().getThreeLetterCode(), chainID).trim() + "\n";
			}
		}
		return out;
	}

	/**
	 * 
	 * @return all backbone atoms as Atom[] not existing atoms are written as
	 *         null and have to be checked in all following methods
	 */
	public Atom[] getBackboneAtoms() {
		AminoAcid tmp;
		Atom[] backbone = new Atom[aminoAcids.length * 4];
		Atom[] temp;
		for (int i = 0; i != aminoAcids.length; i++) {
			tmp = aminoAcids[i];
			temp = tmp.getBackboneAtoms();
			for (int j = 0; j != 4; j++) {
				backbone[i * 4 + j] = temp[j];
			}
		}
		return backbone;
	}

	/**
	 * returns the sequence as a String
	 * @return the amino acid sequence as a String
	 */
	public String getSequence() {
		StringBuilder seq = new StringBuilder();
		for(int i = 0; i < length; i++) {
			seq.append(aminoAcids[i].getName());
		}
		return seq.toString();
	}
}
