package db.mysql;

import bioinfo.proteins.AminoAcid;
import bioinfo.proteins.AminoAcidName;
import bioinfo.proteins.PDBEntry;

public class ImportPDBs {

	public static void main(String[] args) {
		MysqlDBConnection connection = new MysqlDBConnection();
		PDBConnector pdbconnector = new PDBConnector(connection);
		AtomConnector atomconnector = new AtomConnector(connection);
		AAConnector aaconnector = new AAConnector(connection);
		
		AminoAcid[] aminos = {new AminoAcid(AminoAcidName.A),new AminoAcid(AminoAcidName.C)};		
		PDBEntry pdbentry = new PDBEntry("11asB00",aminos);
		pdbconnector.addEntry(pdbentry);
	}

}
