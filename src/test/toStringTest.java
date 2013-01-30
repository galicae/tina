package test;

import java.util.LinkedList;
import java.util.List;

import bioinfo.proteins.PDBEntry;
import bioinfo.proteins.PDBFileReader;

public class toStringTest {

	public static void main(String[] args) {
		PDBFileReader reader = new PDBFileReader("./temp");
		List<PDBEntry> list = new LinkedList<PDBEntry>();
		list = reader.readPdbFolder();
		
		PDBEntry aj = list.get(0);
		int startIndex = 0;
		for (int i = 0; i < 5; i++) {
			System.out.println(aj.getAminoAcid(i).toPDBLineString(startIndex, 'A'));
			startIndex += aj.getAminoAcid(i).getAtomNumber();
		}
	}
}
