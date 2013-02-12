package test;

import bioinfo.proteins.PDBEntry;
import bioinfo.proteins.PDBFileReader;
import bioinfo.proteins.ReProf;
import bioinfo.proteins.SecStructThree;

public class ReprofTest {

	public static void main(String[] args) throws Exception {
		PDBFileReader r = new PDBFileReader("./proteins2");
		r.initSequentialFolderRead();
		PDBEntry pdb = r.nextPdb();

		SecStructThree[] d = ReProf.predictSecStruct(pdb, "/home/galicae/Desktop/", "lib/reprof");
		
		for(int i = 0; i < d.length; i++) {
			System.out.print(d[i]);
		}
	}
}
