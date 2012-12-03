package test;

import bioinfo.proteins.PDBEntry;
import bioinfo.proteins.PDBFileReader;
import util.JMolView;

public class JmolTest {

	public static void main(String[] args) {
		PDBFileReader reader0 = new PDBFileReader();
		PDBEntry entry0 = reader0.readPDBFromFile("/Users/andreseitz/Documents/uni/CIP/gobi/STRUCTURES/1j2xA00.pdb");
		
		PDBFileReader reader1 = new PDBFileReader();
		PDBEntry entry1 = reader1.readPDBFromFile("/Users/andreseitz/Documents/uni/CIP/gobi/STRUCTURES/1wq2B00.pdb");
		
		PDBFileReader reader2 = new PDBFileReader();
		PDBEntry entry2 = reader2.readPDBFromFile("/Users/andreseitz/Documents/uni/CIP/gobi/STRUCTURES/2aduA02.pdb");
		
		
		JMolView jmol = new JMolView();
		jmol.addPDBEntry(entry0,"red");
		jmol.addPDBEntry(entry1,"blue");
		jmol.addPDBEntry(entry2,"green");
	}
}
