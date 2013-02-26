package bioinfo.proteins.corecluster;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileReader;
import java.io.IOException;
import java.io.ObjectInputStream;
import java.io.PrintWriter;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;

import bioinfo.energy.potential.voronoi.VoroPPWrap;
import bioinfo.energy.potential.voronoi.VoronoiData;
import bioinfo.proteins.PDBEntry;
import bioinfo.proteins.PDBFileReader;
import bioinfo.proteins.fragm3nt.ProteinFragment;
import bioinfo.superpos.PDBReduce;

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
		LinkedList<String> ids = new LinkedList<String>();
		readAllIds(ids);
		LinkedList<MultiCurve> allCurves = new LinkedList<MultiCurve>();

		calculateAllMCs(ids, allCurves);

	}

	/**
	 * function to read serialized VectorAnnotation objects
	 * 
	 * @param fileName
	 * @return
	 * @throws ClassNotFoundException
	 * @throws IOException
	 */
	private static VectorAnnotation getVectorAnnotation(String fileName)
			throws ClassNotFoundException, IOException {
		ObjectInputStream objectInputStream = new ObjectInputStream(
				new FileInputStream(new File(fileName)));
		VectorAnnotation result = (VectorAnnotation) objectInputStream
				.readObject();
		objectInputStream.close();
		return result;
	}

	/**
	 * magic function, calculates all pairs of curves in contact
	 * 
	 * @param ids
	 *            the cath ids of the proteins
	 * @param allCurves
	 *            a list to store MultiCurve objects in
	 * @throws Exception
	 */
	private static void calculateAllMCs(LinkedList<String> ids,
			LinkedList<MultiCurve> allCurves) throws IOException, ClassNotFoundException {
		
		MultiCurveEvaluationTree tree = new MultiCurveEvaluationTree();
		
		PDBFileReader re = new PDBFileReader(pdbDir);
		VoroPPWrap voro = new VoroPPWrap(voroBin);
		PrintWriter wr = new PrintWriter("bla");
		for (int i = 0; i < ids.size(); i++) {
			try {
				String currentId = ids.get(i);
				VectorAnnotation curAnnotation = getVectorAnnotation("./vectors/"
						+ currentId + ".vct");
				List<Curve> curves = curAnnotation.getAllVectors();

				PDBEntry e = re.readPDBFromFile(pdbDir + currentId + ".pdb");
				ProteinFragment ef = new ProteinFragment(e.getID(),
						e.getSequenceAsString(), PDBReduce.reduceSinglePDB(e),
						e.length());
				VoronoiData data = new VoronoiData(e.getID());
				data.reducePDB(e);
				voro.decomposite(data);
				HashMap<Integer, HashMap<Integer, Double>> faces = data
						.getFaces();
				for (int j = 0; j < curves.size() - 1; j++) {
					if(curves.get(j).getType() != 'E')
						continue;
//					if(curves.get(j).getSeqLength() < 6)
//						continue;
					for (int k = j + 1; k < curves.size(); k++) {
//						if(k == j)
//							continue;
						if(curves.get(k).getType() != 'E')
							continue;
//						if(curves.get(k).getSeqLength() < 6)
//							continue;
						double contact = contactArea(curves.get(j),
								curves.get(k), faces);
						if (contact > 0) {
							// they are in contact and can be saved!
							MultiCurve mc = new MultiCurve(e.getID(), contact);
							int start = curves.get(j).getPosition().getX();
							int end = curves.get(j).getPosition().getY();
							int[] pos = { start, end };
							mc.addElement(curves.get(j), ef.getPart(pos));
							start = curves.get(k).getPosition().getX();
							end = curves.get(k).getPosition().getY();
							int[] pos2 = { start, end };
							mc.addElement(curves.get(k), ef.getPart(pos2));
							allCurves.add(mc);
							mc.calculateAllAngles();
							String type = curves.get(j).getType() + "" + curves.get(k).getType();
							MultiCurveDataPoint bla = new MultiCurveDataPoint(e.getID(), contact, mc.getTheta(0, 1), type, Math.min(curves.get(j).getSeqLength(), curves.get(k).getSeqLength()));
							tree.insertData(bla);
							wr.printf("%f %f %d\n", mc.getTheta(0, 1) * 180 / Math.PI, contact, Math.min(curves.get(j).getSeqLength(), curves.get(k).getSeqLength()));
						}
					}
				}
			} catch (NullPointerException e) {
				continue;
			}
			tree.printTree();
		}
		wr.close();
	}

	/**
	 * returns the contact area of two curves (vectorized core elements) based
	 * on a voronoi decomposition
	 * 
	 * @param a
	 *            first curve
	 * @param b
	 *            second curve
	 * @param faces
	 *            the distance hashmap produced by the Voronoi decomposition
	 * @return the contact area. Greater than 0 means cores are in contact
	 */
	private static double contactArea(Curve a, Curve b,
			HashMap<Integer, HashMap<Integer, Double>> faces) {
		int norm = 0;
		double totalContactArea = 0;
		for (int i = a.getPosition().getX(); i < a.getPosition().getY(); i++) {
			int j = b.getPosition().getX();
			for (j = b.getPosition().getX(); j < b.getPosition().getY(); j++) {
				try {
					totalContactArea += faces.get(i).get(j);
					norm++;
				} catch (Exception e) {
					continue;
				}
			}
		}
		
		norm = Math.min(a.getSeqLength(), b.getSeqLength());

		return totalContactArea / (norm * 1.0);
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
