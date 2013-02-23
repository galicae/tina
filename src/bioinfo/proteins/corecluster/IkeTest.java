package bioinfo.proteins.corecluster;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileReader;
import java.io.IOException;
import java.io.ObjectInputStream;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;

import bioinfo.energy.potential.voronoi.VoroPPWrap;
import bioinfo.energy.potential.voronoi.VoronoiData;
import bioinfo.proteins.PDBEntry;
import bioinfo.proteins.PDBFileReader;

/**
 * testing the classes necessary for the loop experiment
 * 
 * @author galicae
 * 
 */
public class IkeTest {
	private static String pdbDir = "/home/galicae/Desktop/STRUCTURES/";
	static String voroBin = "./tools/voro++_ubuntuquantal";

	public static void main(String[] args) throws Exception {
		String dsspFile = "all.dssp";
		PDBFileReader re = new PDBFileReader(pdbDir);
		LinkedList<String> ids = new LinkedList<String>();
		VoroPPWrap voro = new VoroPPWrap(voroBin);

		readAllIds(ids);
		LinkedList<MultiCurve> allCurves = new LinkedList<MultiCurve>();

		for (int i = 0; i < ids.size(); i++) {
			String currentId = ids.get(i);
			VectorAnnotation curAnnotation = getVectorAnnotation("./vectors/" + currentId + ".vct");
			List<Curve> curves = curAnnotation.getAllVectors();

			PDBEntry e = re.readPDBFromFile(pdbDir + currentId + ".pdb");
			VoronoiData data = new VoronoiData(e.getID());
			data.reducePDB(e);
			voro.decomposite(data);
			HashMap<Integer, HashMap<Integer, Double>> faces = data.getFaces();
			for (int j = 0; j < curves.size() - 1; j++) {
				for (int k = j + 1; k < curves.size(); k++) {
					double contact = contactArea(curves.get(j), curves.get(k),
							faces);
					if (contact > 0) {
						// they are in contact and can be saved!
						MultiCurve mc = new MultiCurve(e.getID(), contact);
						mc.addElement(curves.get(j));
						mc.addElement(curves.get(k));
						allCurves.add(mc);
						System.out.println(contact + " " + mc.calculateAngle(0, 1));
					}
				}
			}
		}
	}

	private static VectorAnnotation getVectorAnnotation(String fileName) throws ClassNotFoundException, IOException {
		ObjectInputStream objectInputStream = new ObjectInputStream(new FileInputStream(new File(fileName)));
		VectorAnnotation result = (VectorAnnotation) objectInputStream.readObject();
		return result;
	}
	
	/**
	 * 
	 * @param a
	 * @param b
	 * @param faces
	 * @return
	 */
	private static double contactArea(Curve a, Curve b,
			HashMap<Integer, HashMap<Integer, Double>> faces) {
		double totalContactArea = 0;
		for (int i = a.getPosition().getX(); i < a.getPosition().getY(); i++) {
			int j = b.getPosition().getX();
			for (j = b.getPosition().getX(); j < b.getPosition().getY(); j++) {
				try {
					totalContactArea += faces.get(i).get(j);
				} catch (Exception e) {
					continue;
				}
			}
		}

		return totalContactArea;
	}

	/**
	 * shortcut for reading the ids of all pairs in cathscop.ids
	 * 
	 * @param ids
	 */
	private static void readAllIds(LinkedList<String> ids) {
		try {
			BufferedReader r = new BufferedReader(
					new FileReader("cathscop.ids"));
			String line = "";
			while ((line = r.readLine()) != null) {
				ids.add(line.split("\\s+")[0]);
			}
			r.close();
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}
