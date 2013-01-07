package test;

import java.util.LinkedList;
import java.util.List;

import bioinfo.proteins.PDBEntry;
import bioinfo.proteins.PDBFileReader;
import bioinfo.proteins.fragm3nt.Fragmenter;

public class Fr4gmentTest {
	
	public static void main(String[] args) {
		PDBFileReader reader = new PDBFileReader("./proteins/");
		
		List<PDBEntry> files = new LinkedList<PDBEntry>();
		LinkedList<char[]> pList = new LinkedList<char[]>();
		PDBEntry pdb1 = reader.readPDBFromFile(args[0]);
		
		files.add(pdb1);
		System.out.println("starting disassembly");
		for(PDBEntry e: files) {
			Fragmenter.disassemble(e, pList, 7);
		}
		System.out.println("starting alignment");
		
	}

}
