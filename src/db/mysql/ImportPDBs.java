package db.mysql;

import bioinfo.pdb.PDBFile;
import bioinfo.proteins.AminoAcid;
import bioinfo.proteins.AminoAcidName;
import bioinfo.proteins.PDBEntry;
import bioinfo.proteins.PDBFileReader;

public class ImportPDBs {

	public static void main(String[] args) {
		MysqlDBConnection connection = new MysqlDBConnection();
		PDBConnector pdbconnector = new PDBConnector(connection);
		AtomConnector atomconnector = new AtomConnector(connection);
		AAConnector aaconnector = new AAConnector(connection);
		
		PDBFileReader pdbreader = new PDBFileReader(args[0]);
		PDBFile.downloadPDB("2ADU", "./");
		PDBEntry test = pdbreader.readFromFolderById("1J2xB00");
		for (int i = 0; i < test.length(); i++) {
			System.out.println(test.getAminoAcid(i).getResIndex());
		}
		
		
		
//		AminoAcid[] aminos = {new AminoAcid(AminoAcidName.A),new AminoAcid(AminoAcidName.C)};		
//		PDBEntry pdbentry = new PDBEntry("11asB00",aminos);
//		pdbconnector.addEntry(pdbentry);
	}

}
