package bioinfo.superpos;

import bioinfo.alignment.Alignment;
import bioinfo.alignment.SequenceAlignmentFileReader;
import bioinfo.proteins.PDBEntry;
import bioinfo.proteins.PDBFileReader;

/**
 * this is a naive take to exercise 2 of assignment 1
 * 
 * @author nikos
 * 
 */
public class TMCollective {
	
	public TMCollective() {

	}

	public PDBEntry[] createTMInput(String align, String pFile, String qFile) {
		SequenceAlignmentFileReader ali = new SequenceAlignmentFileReader(align);
		Alignment alignment = ali.readAlignments().get(0);
		PDBFileReader reader = new PDBFileReader();

		//ArrayList<AminoAcid> pAminoList = new ArrayList<AminoAcid>();
		//ArrayList<AminoAcid> qAminoList = new ArrayList<AminoAcid>();

		PDBEntry p = reader.readPDBFromFile(pFile);
		PDBEntry q = reader.readPDBFromFile(qFile);

		// return only CA atoms of amino acids in PDB form
		PDBEntry[] result = PDBReduce.reducePDB(alignment, p, q);
		return result;
	}

	public PDBEntry[] createTMInput(Alignment alignment, PDBEntry pFile,
			PDBEntry qFile) {

		// return only CA atoms of amino acids in PDB form
		PDBEntry[] result = PDBReduce.reducePDB(alignment, pFile, qFile);
		return result;
	}
}
