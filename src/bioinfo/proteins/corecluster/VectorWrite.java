package bioinfo.proteins.corecluster;

import java.io.BufferedReader;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.IOException;
import java.io.ObjectOutputStream;
import java.util.LinkedList;
import java.util.Map;


public class VectorWrite {
	private static String pdbDir = "/home/galicae/Desktop/STRUCTURES/";
	static String voroBin = "./tools/voro++_ubuntuquantal";

	public static void main(String[] args) throws Exception {
		String dsspFile = "all.dssp";
		LinkedList<String> ids = new LinkedList<String>();

		readAllIds(ids);

		VectorAnalyzer vv = new VectorAnalyzer(pdbDir, dsspFile);
		Map<String, VectorAnnotation> refined = vv.letThereBeDragons();

		for (int i = 0; i < ids.size(); i++) {
			String currentId = ids.get(i);
			VectorAnnotation curAnnotation = refined.get(currentId);
			writeVectorAnnotation(curAnnotation, "./vectors/");
		}
	}

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
	
	private static void writeVectorAnnotation(VectorAnnotation vectorAnnotation, String outputDir) throws IOException {
		ObjectOutputStream objectOutputStream = new ObjectOutputStream(new FileOutputStream(outputDir + vectorAnnotation.getId() + ".vct"));
		objectOutputStream.writeObject(vectorAnnotation);
		objectOutputStream.close();
	}
}
