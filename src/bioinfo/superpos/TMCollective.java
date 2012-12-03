package bioinfo.superpos;

import java.util.ArrayList;

import test.tempAlignmentReader;
import bioinfo.alignment.Alignment;
import bioinfo.proteins.AminoAcid;
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
		tempAlignmentReader ali = new tempAlignmentReader(align);
		Alignment alignment = ali.readAlignment();
		PDBFileReader reader = new PDBFileReader();

		//ArrayList<AminoAcid> pAminoList = new ArrayList<AminoAcid>();
		//ArrayList<AminoAcid> qAminoList = new ArrayList<AminoAcid>();

		PDBEntry p = reader.readPDBFromFile(pFile);
		PDBEntry q = reader.readPDBFromFile(qFile);

		// return only CA atoms of amino acids in PDB form
		PDBReduce.reducePDB(alignment, p, q);
		PDBEntry[] result = new PDBEntry[2];
		result[0] = p;
		result[1] = q;
		return result;
	}
}
