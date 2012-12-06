package test;

import bioinfo.alignment.SequenceAlignment;
import bioinfo.alignment.SequenceAlignmentFileReader;
import bioinfo.proteins.PDBEntry;
import bioinfo.proteins.PDBFileReader;
import bioinfo.superpos.Kabsch;
import bioinfo.superpos.PDBReduce;
import bioinfo.superpos.TMMain;
import bioinfo.superpos.Transformation;

public class TMPipeline {
	public static void main(String[] args) throws Exception {

		PDBFileReader reader = new PDBFileReader();
		PDBEntry pdb1 = reader
				.readPDBFromFile("C:/Users/nikos/Desktop/STRUCTURES/1muzA00.pdb");
		PDBEntry pdb2 = reader
				.readPDBFromFile("C:/Users/nikos/Desktop/STRUCTURES/1k4uS00.pdb");

//		System.out.println(pdb1.getAtomSectionAsString());
		SequenceAlignmentFileReader aliReader = new SequenceAlignmentFileReader(
				args[2]);
		SequenceAlignment alignment = aliReader.readAlignments().get(0);

		TMMain main = new TMMain();
//		Transformation trOr = main.calculateTransformation(args[2],"1muzA00.pdb","1k4uS00.pdb");
		Transformation tr1 = main.calculateTransformation(alignment, pdb1, pdb2);
		System.out.println(tr1.getRmsd());
//		System.out.println(trOr.getTmscore());
	}
}
