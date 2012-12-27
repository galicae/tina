package bioinfo.proteins.framg3nt;

import java.util.LinkedList;

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
	public static LinkedList<ProteinFragment> crunch(PDBEntry pdb,
			LinkedList<ProteinFragment> l, int fLength) {
		for (int i = 0; i < pdb.length() - fLength; i++) {
			double[][] temp = new double[4 * fLength][3];
			for (int j = i; j < i + fLength; j++) {				
				temp[(j-i) * 4 + 0] = pdb.getAminoAcid(j).getAtomByType(AtomType.N)
						.getPosition();
				temp[(j-i) * 4 + 1] = pdb.getAminoAcid(j)
						.getAtomByType(AtomType.CA).getPosition();
				temp[(j-i) * 4 + 2] = pdb.getAminoAcid(j).getAtomByType(AtomType.O)
						.getPosition();
				temp[(j-i) * 4 + 3] = pdb.getAminoAcid(j).getAtomByType(AtomType.C)
						.getPosition();
			}
			ProteinFragment tempFrag = new ProteinFragment(pdb.getID() + "_"
					+ i, temp, i, fLength);
			l.add(tempFrag);
		}
		return l;
	}
}
