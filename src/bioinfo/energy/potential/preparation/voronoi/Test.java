package bioinfo.energy.potential.preparation.voronoi;

import java.io.IOException;

import bioinfo.proteins.PDBEntry;
import bioinfo.proteins.PDBFileReader;

public class Test {

	/**
	 * @param args
	 */
	public static void main(String[] args) {

		PDBFileReader reader = new PDBFileReader();
		PDBEntry pdb = reader.readPDBFromFile("/Users/andreseitz/Documents/uni/CIP/gobi/STRUCTURES/1j2xA00.pdb");
		VoroPrepare vprep = new VoroPrepare();
		VoroPrepType type = VoroPrepType.CC;
		VoronoiData data = vprep.reducePDB(type, pdb);
		//VoroPPWrap voro = new VoroPPWrap();
		//data = voro.decomposite(data);
		
		
	}

}
