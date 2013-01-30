package db.mysql;

//import bioinfo.pdb.PDBFile;
import bioinfo.proteins.AtomType;
import bioinfo.proteins.PDBEntry;
import bioinfo.proteins.PDBFileReader;

public class ImportPDBs {

	public static void main(String[] args) {
		LocalConnection connection = new LocalConnection();
		PDBConnector pdbconnector = new PDBConnector(connection);
		
		PDBFileReader pdbreader = new PDBFileReader(args[0]);
		
		String pdb = "1j2xA00";
		
		if(!pdbconnector.pdbExist(pdb)){
			PDBEntry test = pdbreader.readFromFolderById(pdb);
			pdbconnector.addEntry(test);
		} else{
			System.out.println("PDB '"+pdb+"' existiert bereits in DB.");
		}
		
		PDBEntry out = pdbconnector.getPDB(pdb);
		System.out.println(out.getAminoAcid(7).getAtomByType(AtomType.CA).getPosition()[0]);
	}

}
