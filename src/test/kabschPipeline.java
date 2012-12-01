package test;

import bioinfo.alignment.SequenceAlignment;
import bioinfo.proteins.PDBEntry;
import bioinfo.proteins.PDBFileReader;
import bioinfo.superpos.*;

public class kabschPipeline {

	public static void main(String[] args) {
		
		PDBFileReader reader = new PDBFileReader();
		PDBEntry pdb1 = reader.readPDBFromFile(args[1]);
		PDBEntry pdb2 = reader.readPDBFromFile(args[2]);
		tempAlignmentReader aliReader = new tempAlignmentReader(args[0]);
		SequenceAlignment alignment = (SequenceAlignment)aliReader.readAlignment();
		
		
		double[][][] reducedPdbs = PDBReduce.reduce(alignment, pdb1, pdb2);
		Transformation tr = Kabsch.calculateTransformation(reducedPdbs);
		PDBEntry superposedPdb2 = tr.transform(pdb2);
		
		System.out.println(superposedPdb2.getID());
		System.out.println(tr.getRmsd());
	}
}
