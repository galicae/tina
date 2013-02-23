package test;

import java.util.HashMap;

import bioinfo.energy.potential.voronoi.VoroPPWrap;
import bioinfo.energy.potential.voronoi.VoronoiData;
import bioinfo.proteins.PDBEntry;
import bioinfo.proteins.PDBFileReader;

public class VoroTest {
	public static void main(String[] args) {
		PDBFileReader r = new PDBFileReader("/home/galicae/Desktop/STRUCTURES/");
		PDBEntry e = r.readFromFolderById("1j2xA00");
		
		String voroBin = "./tools/voro++_ubuntuquantal";
		
		VoroPPWrap voro = new VoroPPWrap(voroBin);
		VoronoiData data = new VoronoiData(e.getID());
		data.reducePDB(e);
		voro.decomposite(data);
		
		HashMap<Integer, HashMap<Integer, Double>> faces = data.getFaces();
	}
}
