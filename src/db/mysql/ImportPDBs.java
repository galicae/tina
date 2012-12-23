package db.mysql;

import bioinfo.pdb.PDBFile;
import bioinfo.proteins.AtomType;
import bioinfo.proteins.PDBEntry;
import bioinfo.proteins.PDBFileReader;

public class ImportPDBs {

	public static void main(String[] args) {
		MysqlDBConnection connection = new MysqlDBConnection();
		PDBConnector pdbconnector = new PDBConnector(connection);
		
		PDBFileReader pdbreader = new PDBFileReader(args[0]);
		//PDBFile.downloadPDB("2ADU", "./");
		PDBEntry test = pdbreader.readFromFolderById("1J2xB00");		
				
		PDBEntry out = pdbconnector.getPDB("1J2X");
		System.out.println(out.getAminoAcid(7).getAtomByType(AtomType.CA).getPosition()[0]);
	}

}
