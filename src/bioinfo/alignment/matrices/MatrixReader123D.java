package bioinfo.alignment.matrices;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;

import util.Bio;

/**
 * Reader class for all matrices that 123D needs apart from subsitution matrix
 * 
 * @author gobi_4
 * @date December 5, 2012
 * 
 */
public class MatrixReader123D {

	/**
	 * reads a file containing the secondary structure preferences for all amino
	 * acids given in 3-letter code, with a caps header (ALPHA-BETA-OTHER) over
	 * each section and an UNK entry
	 * 
	 * @param filename
	 * @return the secondary structure preference as a double array
	 */
	public double[][] readSecStructPref(String filename) {
		double[][] matrix = new double[3][27];
		try {
			BufferedReader in = new BufferedReader(new FileReader(filename));
			String line;
			String[] content = new String[2];
			while ((line = in.readLine()) != null) {
				if (line.contains("ALPHA")) {
					while (!(line = in.readLine()).startsWith("UNK")
							&& (line = in.readLine()) != null) {
						content = line.split("\\s+");
						char c = Bio.codeTranslate(content[0]);
						matrix[0][c - 65] = Double.parseDouble(content[1]);
					}
					matrix[0][26] = 0;
				} else if (line.contains("BETA")) {
					while (!(line = in.readLine()).startsWith("UNK")
							&& (line = in.readLine()) != null) {
						content = line.split("\\s+");
						char c = Bio.codeTranslate(content[0]);
						matrix[1][c - 65] = Double.parseDouble(content[1]);
					}
					matrix[1][26] = 0;
				} else if (line.contains("OTHER")) {
					while (!(line = in.readLine()).startsWith("UNK")
							&& (line = in.readLine()) != null) {
						content = line.split("\\s+");
						char c = Bio.codeTranslate(content[0]);
						matrix[2][c - 65] = Double.parseDouble(content[1]);
					}
					matrix[2][26] = 0;
				}
			}
			in.close();
		} catch (IOException e) {
			System.out.println("No Input (secondary structure preference)!");
		}
		return matrix;
	}

}
