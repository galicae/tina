package bioinfo.proteins.fr4gment;

import bioinfo.proteins.PDBEntry;
import bioinfo.proteins.PDBFileReader;
import bioinfo.superpos.PDBReduce;

public class ModellerTest {
	// /home/proj/biosoft/PROTEINS/scripts/model_full.sh -i
	// /home/p/papadopoulos/Desktop/test.ali -S
	// /home/proj/biosoft/PROTEINS/CATHSCOP/STRUCTURES/ -o 1j2xA00.pdb

	public static void main(String[] args) {
		PDBFileReader read = new PDBFileReader();
		PDBEntry model = read
				.readPDBFromFile("/home/p/papadopoulos/modeller/1xpxA00.pdb");

		double[][] modelCoordinates = PDBReduce.reduceSinglePDB(model);
		System.out.println(modelCoordinates);
	}

}
