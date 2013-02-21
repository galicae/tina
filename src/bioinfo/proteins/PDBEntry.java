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

import bioinfo.Sequence;
import bioinfo.StructuralData;
import bioinfo.alignment.Alignable;

/**
 * @author gobi_4
 * @lastchange 2013-02-18
 */
public class PDBEntry implements Alignable, Serializable, StructuralData {

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
	 * stanard constructor for PDBEntry XXXXA00
	 * 
	 * @param id
	 *            the ID of the PDB Entry (XXXX)
	 * @param chainID
	 *            the chain ID (A)
	 * @param chainIDNum
	 *            the number of the chain (00)
	 * @param aminoAcids
	 *            Array of AminoAcids
	 */
	public PDBEntry(String id, char chainID, int chainIDNum,
			AminoAcid[] aminoAcids) {
		this.id = id;
		this.chainID = chainID;
		this.chainIDNum = chainIDNum;
		this.setAminoAcids(aminoAcids);
	}

	/**
	 * Constructor for an a PDBEntry with an Array of AminoAcids
	 * 
	 * @param id
	 *            the id of the PDBEntry (XXXX or XXXXA00)
	 * @param aminoAcids
	 *            Array of AminoAcids
	 */
	public PDBEntry(String id, AminoAcid[] aminoAcids) {
		if (id.length() == 4) {
			this.id = id;
			this.chainID = 'A';
			this.chainIDNum = 0;
		} else {
			this.id = id.substring(0, 4);
			this.chainID = id.charAt(4);
			this.chainIDNum = Integer.valueOf(id.substring(5, 7));
		}
		this.setAminoAcids(aminoAcids);
	}

	/**
	 * Constructor for an the PDBEntry without any AminoAcids
	 * 
	 * @param id
	 *            the id of the PDBEntry (XXXX or XXXXA00)
	 */
	public PDBEntry(String id) {
		this(id, new AminoAcid[0]);
	}

	/**
	 * Constructor for an the PDBEntry with a List of AminoAcids
	 * 
	 * @param id
	 *            the id of the PDBEntry (XXXX or XXXXA00)
	 * @param aminoAcidList
	 *            List of AminoAcids
	 */
	public PDBEntry(String id, List<AminoAcid> aminoAcidList) {
		this(id, aminoAcidList.toArray(new AminoAcid[0]));
	}

	/**
	 * Constructor for an an PDBEntry. don't know what this is good for. Doesn't
	 * even seem to be needed.
	 * 
	 * @param id
	 *            the id of the PDBEntry (XXXX or XXXXA00)
	 * @param aminoAcidList
	 *            List of AminoAcids
	 * @param real
	 *            ?????
	 */
	public PDBEntry(String id, List<AminoAcid> aminoAcidList, boolean real) {
		if (id.length() == 4) {
			this.id = id;
			this.chainID = 'A';
			this.chainIDNum = 0;
		} else {
			this.id = id.substring(0, 4);
			this.chainID = id.charAt(4);
			this.chainIDNum = Integer.valueOf(id.substring(5));
			this.setAminoAcids(aminoAcidList.toArray(new AminoAcid[0]));
		}
	}

	/**
	 * Constructor for a new PDBEntry which has the same Values as the given one
	 * 
	 * @param arg
	 *            PDBEntry with the Parameters to copy from
	 */
	public PDBEntry(PDBEntry arg) {
		this(arg.id, arg.chainID, arg.chainIDNum, arg.aminoAcids);
	}

	/**
	 * return complete id (xxxxA00)
	 * @return
	 */
	public String getID(){
		return id+chainID+getChainIDNumAsString();
	}
	
	/**
	 * @return the ID without chain and num of the PDBEntry (xxxx)
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
			result[i] = getAminoAcids()[i];
		}
		result[length] = arg;
		setAminoAcids(result);
		length++;
	}

	/**
	 * @param n
	 * @return the n-th AminoAcid
	 * @throws NullPointerException
	 *             if n is out of bounds
	 */
	public AminoAcid getAminoAcid(int n) {
		if (n >= aminoAcids.length) {
			throw new IndexOutOfBoundsException(
					"This AminoAcid doesn't exist in this PDBEntry!");
		}
		return getAminoAcids()[n];
	}

	/**
	 * @return number of AminoAcids in this PDBEntry
	 */
	@Override
	public int length() {
		return getAminoAcids().length;
	}

	/**
	 * @param n
	 *            integer defining position of wanted AminoAcid in PDB Entry
	 * @return the n-th AminoAcid
	 */
	@Override
	public AminoAcid getComp(int n) {
		if (n >= aminoAcids.length) {
			throw new IndexOutOfBoundsException(
					"This AminoAcid doesn't exist in this PDBEntry!");
		}
		return getAminoAcids()[n];
	}

	/**
	 * returns the AmioAcids as a (human readable) String like in a PDB file
	 * 
	 * @return all AminoAcids as a (human readable) String
	 */
	public String getAtomSectionAsString() {
		int lineCounter = 0;
		String out = "";
		for (int i = 0; i != getAminoAcids().length; i++) {
			for (int j = 0; j != getAminoAcids()[i].getAtomNumber(); j++) {
				// TODO Some Entries have no 7-letter ID (1TIMA00) but four
				// letters (1TIM)
				// ==> make some Method for IDs without chainID!
				out += getAminoAcids()[i]
						.getAtom(j)
						.toString(
								lineCounter++,
								i,
								getAminoAcids()[i].getName()
										.getThreeLetterCode(), chainID).trim()
						+ "\n";
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
		Atom[] backbone = new Atom[getAminoAcids().length * 4];
		Atom[] temp;
		for (int i = 0; i != getAminoAcids().length; i++) {
			tmp = getAminoAcids()[i];
			temp = tmp.getBackboneAtoms();
			for (int j = 0; j != 4; j++) {
				backbone[i * 4 + j] = temp[j];
			}
		}
		return backbone;
	}

	/**
	 * returns the sequence as a String
	 * 
	 * @return the amino acid sequence as a String
	 */
	public String getSequenceAsString() {
		StringBuilder seq = new StringBuilder();
		for (int i = 0; i < length; i++) {
			seq.append(getAminoAcids()[i].getName());
		}
		return seq.toString();
	}

	/**
	 * @return the aminoAcids
	 */
	public AminoAcid[] getAminoAcids() {
		return aminoAcids;
	}

	/**
	 * sets aminoAcid and updates this.length
	 * 
	 * @param aminoAcids
	 *            the aminoAcids to set
	 */
	public void setAminoAcids(AminoAcid[] aminoAcids) {
		this.aminoAcids = aminoAcids;
		this.length = this.aminoAcids.length;
	}

	/**
	 * 
	 * @return the Sequence that is depicted by this PDBEntry
	 */
	public Sequence getSequence() {
		return new Sequence(getID(), getSequenceAsString());
	}

	/**
	 * 
	 * @returnchainIDNum as two-lengthed String
	 */
	private String getChainIDNumAsString() {
		String result = Integer.toString(chainIDNum);
		if (chainIDNum < 10) {
			result = "0" + result;
		}
		return result;
	}
	
	public String toString() {
		return "PDBEntry: "+this.getSequenceAsString();
	}

}
