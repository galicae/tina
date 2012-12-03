package test;

import bioinfo.proteins.PDBEntry;
import bioinfo.proteins.PDBFileReader;

public class toStringTest {
	
	public static void main(String[] args) {
		PDBFileReader re = new PDBFileReader();
		PDBEntry ne = re.readPDBFromFile(args[0]);
		
		System.out.println(ne.getAminoAcid(0).getAtom(0).toString());
	}
}
