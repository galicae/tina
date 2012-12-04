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
		PDBEntry pdb1 = reader.readPDBFromFile(args[0]);
		PDBEntry pdb2 = reader.readPDBFromFile(args[1]);
		PDBEntry pdb2Copy = reader.readPDBFromFile(args[1]);
		SequenceAlignmentFileReader aliReader = new SequenceAlignmentFileReader(args[2]);
		SequenceAlignment alignment = aliReader.readAlignments().get(0);
		
		TMMain main = new TMMain();
		Transformation tr = main.calculateTransformation(alignment, pdb1, pdb2);
		double[][][] reducedPdbs = PDBReduce.reduce(alignment, pdb1, pdb2);
		
		PDBEntry superposedPdb2 = tr.transform(pdb2);
		
		System.out.println(superposedPdb2.getID());
		System.out.println(tr.getRmsd());
		
		tr = Kabsch.calculateTransformation(reducedPdbs);
		superposedPdb2 = tr.transform(pdb2Copy);
		
		System.out.println(superposedPdb2.getID());
		System.out.println(tr.getRmsd());
	}
}
