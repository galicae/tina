package db.mysql;

//import bioinfo.pdb.PDBFile;
import java.util.HashMap;
import java.util.Map.Entry;

import bioinfo.alignment.kerbsch.temp.SeqLibrary;
import bioinfo.proteins.DSSPEntry;
import bioinfo.proteins.DSSPFileReader;
import bioinfo.proteins.PDBEntry;
import bioinfo.proteins.PDBFileReader;

public class ImportPDBs {

	public static void main(String[] args) {
		LocalConnection connection = new LocalConnection();
		DBConnector pdbconnector = new DBConnector(connection);
		
		HashMap<String,char[]> seqlib = SeqLibrary.read("../GoBi_old/full_domains.seqlib");
		PDBFileReader pdbreader = new PDBFileReader(args[0]);
		DSSPFileReader dsspreader = new DSSPFileReader(args[1]);
		PDBEntry pdbtest;
		DSSPEntry dssptest;
		
		for(Entry<String,char[]> e : seqlib.entrySet()){
			pdbtest = pdbreader.readFromFolderById(e.getKey());
			pdbconnector.addPDBEntry(pdbtest);
			dssptest = dsspreader.readFromFolderById(e.getKey());
			pdbconnector.addDSSPEntry(dssptest);
		}
		
//		String pdb = "1bus000";
//		String pdb2 = "1e3hA02";
	
		
//		pdbtest = pdbreader.readFromFolderById(pdb);
//		pdbconnector.addPDBEntry(pdbtest);	

//		dssptest = dsspreader.readFromFolderById(pdb2);
//		pdbconnector.addDSSPEntry(dssptest);
		
//		DSSPEntry dssp = pdbconnector.getDSSP(pdb);
//		System.out.println(dssp.getLength());
	}

}
