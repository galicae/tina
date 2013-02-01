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
		
		long start;
		long end;
		long steadyDB = 0;
		long steadyFile = 0;
		int count = 0;
		
//		for(Entry<String,char[]> e : seqlib.entrySet()){
//			pdbtest = pdbreader.readFromFolderById(e.getKey());
//			pdbconnector.addPDBEntry(pdbtest);
//			dssptest = dsspreader.readFromFolderById(e.getKey());
//			pdbconnector.addDSSPEntry(dssptest);
//			
//			start = System.currentTimeMillis();
//			dssptest = pdbconnector.getDSSP(e.getKey());
//			end = System.currentTimeMillis();
//			steadyDB += (end-start);
//			
//			start = System.currentTimeMillis();
//			dssptest = dsspreader.readFromFolderById(e.getKey());
//			end = System.currentTimeMillis();
//			steadyFile += (end-start);
//			
//			System.out.println(++count);
//		}
//		System.out.println("time for file (DSSP): "+ (steadyFile/9303));
//		System.out.println("time for DB (DSSP): "+ (steadyDB/9303));
		
		
		String pdb = "2f0aD00";
//		String pdb2 = "1e3hA02";
	


		
//		start = System.currentTimeMillis();
//		dssptest = dsspreader.readFromFolderById(pdb);
//		end = System.currentTimeMillis();
//		
//		System.out.println("time for file (DSSP): "+ (end-start));
//		
//		start = System.currentTimeMillis();
//		dssptest = pdbconnector.getDSSP(pdb);
//		end = System.currentTimeMillis();
//		
//		System.out.println("time for DB (DSSP): "+ (end-start));
	
		//		pdbconnector.addPDBEntry(pdbtest);	

//		dssptest = dsspreader.readFromFolderById(pdb2);
//		pdbconnector.addDSSPEntry(dssptest);
		
		PDBEntry pdbentry = pdbconnector.getPDB(pdb);
		DSSPEntry dssp = pdbconnector.getDSSP(pdb);
		System.out.println(dssp.getLength());
		System.out.println(pdbentry.length());
		System.out.println();
	}

}
