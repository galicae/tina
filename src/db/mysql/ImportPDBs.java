package db.mysql;

import bioinfo.pdb.PDBFile;
import bioinfo.proteins.PDBEntry;
import bioinfo.proteins.PDBFileReader;

public class ImportPDBs {

	public static void main(String[] args) {
		MysqlDBConnection connection = new MysqlDBConnection();
		PDBConnector pdbconnector = new PDBConnector(connection);
		
		PDBFileReader pdbreader = new PDBFileReader(args[0]);
		PDBFile.downloadPDB("2ADU", "./");
		PDBEntry test = pdbreader.readFromFolderById("1J2xB00");		
				
		pdbconnector.addEntry(test);
	}

}
