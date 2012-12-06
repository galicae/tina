package test;

import bioinfo.proteins.PDBEntry;
import bioinfo.proteins.PDBFileReader;

public class toStringTest {

	public static void main(String[] args) {
		PDBFileReader re = new PDBFileReader();
		PDBEntry ne = re.readPDBFromFile("C:/Users/nikos/Desktop/STRUCTURES/1a0aA00.pdb");

		char c = 65;
		System.out.println(ne.getAminoAcid(0).getAtom(0).toTMString());
		System.out.println(ne.getAminoAcid(0).getAtom(1).toString(0, 0, "ALA", c));
		System.out.println(ne.getAminoAcid(0).getAtom(2).toTMString());
		System.out.println(ne.getAminoAcid(0).getAtom(3).toTMString());
		System.out.println(ne.getAminoAcid(0).getAtom(4).toTMString());
	}
}
