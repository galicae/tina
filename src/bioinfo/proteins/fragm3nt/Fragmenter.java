package bioinfo.proteins.fragm3nt;

import java.util.LinkedList;
import java.util.List;

import bioinfo.proteins.AminoAcid;
import bioinfo.proteins.AtomType;
import bioinfo.proteins.PDBEntry;

public class Fragmenter {

	/**
	 * goes through a PDB Entry and reads protein fragments using the backbone
	 * atoms
	 * 
	 * @param pdb
	 *            the entry to crunch
	 * @param l
	 *            a list of protein fragments
	 * @param fLength
	 *            the desired length of protein fragments
	 * @return the new list of fragments
	 */
	public static LinkedList<ProteinFragment> crunchBackboneN(PDBEntry pdb,
			LinkedList<ProteinFragment> l, int fLength) {
		try {
			AminoAcid tempAA = new AminoAcid("ALA", 0);
			for (int i = 0; i < pdb.length() - fLength; i++) {
				double[][] temp = new double[4 * fLength][3];
				for (int j = i; j < i + fLength; j++) {
					tempAA = pdb.getAminoAcid(j);
					temp[(j - i) * 4 + 0] = tempAA.getAtomByType(AtomType.N)
							.getPosition();
					temp[(j - i) * 4 + 1] = tempAA.getAtomByType(AtomType.CA)
							.getPosition();
					temp[(j - i) * 4 + 2] = tempAA.getAtomByType(AtomType.O)
							.getPosition();
					temp[(j - i) * 4 + 3] = tempAA.getAtomByType(AtomType.C)
							.getPosition();
				}
				ProteinFragment tempFrag = new ProteinFragment(pdb.getID()
						+ "_" + i, temp, i, fLength);
				l.add(tempFrag);
			}
			return l;
		} catch (Exception e) {
			System.out.println("Entry " + pdb.getID()
					+ " probably has incomplete records.");
		}
		return null;
	}

	/**
	 * alternative version that returns the sequence along with the other stuff.
	 * 
	 * @param pdb
	 *            the entry to crunch
	 * @param l
	 *            a list of protein fragments
	 * @param fLength
	 *            the desired length of protein fragments
	 * @return the new list of fragments
	 */
	public static List<ProteinFragment> crunchBackboneSeq(PDBEntry pdb,
			List<ProteinFragment> l, int fLength) {
		try {
			AminoAcid tempAA = new AminoAcid("ALA", 0);
			String seq = "";
			for (int i = 0; i < pdb.length() - fLength; i++) {
				double[][] temp = new double[4 * fLength][3];
				seq = "";
				for (int j = i; j < i + fLength; j++) {
					tempAA = pdb.getAminoAcid(j);
					seq += tempAA.getName().getOneLetterCode();
					temp[(j - i) * 4 + 0] = tempAA.getAtomByType(AtomType.N)
							.getPosition();
					temp[(j - i) * 4 + 1] = tempAA.getAtomByType(AtomType.CA)
							.getPosition();
					temp[(j - i) * 4 + 2] = tempAA.getAtomByType(AtomType.O)
							.getPosition();
					temp[(j - i) * 4 + 3] = tempAA.getAtomByType(AtomType.C)
							.getPosition();
				}
				ProteinFragment tempFrag = new ProteinFragment(pdb.getID()
						+ "_" + i, seq, temp, i, fLength);
				l.add(tempFrag);
			}
			return l;
		} catch (Exception e) {
			System.out.println("Entry " + pdb.getID()
					+ " probably has incomplete records.");
		}
		return null;
	}

	/**
	 * for whenever we need to only make sequence fragments
	 * 
	 * @param pdb
	 *            the entry to disassemble
	 * @param l
	 *            a list of protein fragments
	 * @param fLength
	 *            the desired length of protein fragments
	 * @return the new list of fragments
	 */
	public static LinkedList<char[]> disassemble(PDBEntry pdb,
			LinkedList<char[]> l, int fragmentLength) {
		char[] pdbChar = new char[pdb.length()];
		for (int i = 0; i < pdb.length(); i++) {
			pdbChar[i] = pdb.getAminoAcid(i).getName().getOneLetterCode()
					.charAt(0);
		}

		for (int i = 0; i < pdb.length() - fragmentLength; i++) {
			char[] temp = new char[fragmentLength];
			for (int j = 0; j < fragmentLength; j++) {
				temp[j] = pdbChar[i + j];
			}
			l.add(temp);
		}
		return l;
	}
}
