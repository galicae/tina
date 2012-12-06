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

	/**
	 * wrapper function that creates the input for a TMscore run from files
	 * 
	 * @param align
	 *            the path to an alignment file
	 * @param pFile
	 *            the path to a PDB file for the first PDB entry
	 * @param qFile
	 *            the path to a PDB file for the second PDB entry
	 * @return two "pruned" PDB entries containing only the Ca atoms of aligned
	 *         amino acids
	 */
	public PDBEntry[] createTMInput(String align, String pFile, String qFile) {
		SequenceAlignmentFileReader ali = new SequenceAlignmentFileReader(align);
		Alignment alignment = ali.readAlignments().get(0);
		PDBFileReader reader = new PDBFileReader();

		// ArrayList<AminoAcid> pAminoList = new ArrayList<AminoAcid>();
		// ArrayList<AminoAcid> qAminoList = new ArrayList<AminoAcid>();

		PDBEntry p = reader.readPDBFromFile(pFile);
		PDBEntry q = reader.readPDBFromFile(qFile);

		// return only CA atoms of amino acids in PDB form
		PDBEntry[] result = PDBReduce.reducePDB(alignment, p, q);
		return result;
	}

	/**
	 * the same as above; this function isn't really needed, but was written
	 * nevertheless in order to preserve the general workflow of the TM pipeline
	 * 
	 * @param alignment an alignment object
	 * @param pFile PDBEntry of the first structure
	 * @param qFile PDBEntry of the second structure
	 * @return 
	 */
	public PDBEntry[] createTMInput(Alignment alignment, PDBEntry pFile,
			PDBEntry qFile) {

		// return only CA atoms of amino acids in PDB form
		PDBEntry[] result = PDBReduce.reducePDB(alignment, pFile, qFile);
		return result;
	}
}
