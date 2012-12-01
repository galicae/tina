package test;

import bioinfo.alignment.SequenceAlignment;
import bioinfo.proteins.PDBEntry;
import bioinfo.proteins.PDBFileReader;
import bioinfo.superpos.*;

public class kabschPipeline {

	public static void main(String[] args) {
		
		PDBFileReader reader = new PDBFileReader();
		PDBEntry pdb1 = reader.readPDBFromFile(args[0]);
		PDBEntry pdb2 = reader.readPDBFromFile(args[1]);
		AlignmentReader aliReader = new AlignmentReader(args[2]);
		SequenceAlignment alignment = aliReader.readAlignment();
		
		
		double[][][] reducedPdbs = PDBReduce.reduce(alignment, pdb1, pdb2);
		Transformation tr = Kabsch.calculateTransformation(reducedPdbs);
		PDBEntry superposedPdb2 = tr.transform(pdb2);
	}
}
