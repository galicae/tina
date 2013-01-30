package db.mysql;

//import bioinfo.pdb.PDBFile;
import bioinfo.proteins.DSSPEntry;
import bioinfo.proteins.DSSPFileReader;
import bioinfo.proteins.PDBEntry;
import bioinfo.proteins.PDBFileReader;

public class ImportPDBs {

	public static void main(String[] args) {
		LocalConnection connection = new LocalConnection();
		PDBConnector pdbconnector = new PDBConnector(connection);
		
		PDBFileReader pdbreader = new PDBFileReader(args[0]);
		DSSPFileReader dsspreader = new DSSPFileReader(args[1]);
		
		String pdb = "1j2xA00";
	
		PDBEntry pdbtest = pdbreader.readFromFolderById(pdb);
		DSSPEntry dssptest = dsspreader.readFromFolderById(pdb);
		//pdbconnector.addPDBEntry(pdbtest);
		//pdbconnector.addDSSPEntry(dssptest);

		DSSPEntry out = pdbconnector.getDSSP(pdb);
		System.out.println(out.getID());
		System.out.println(out.getLength());
	}

}
