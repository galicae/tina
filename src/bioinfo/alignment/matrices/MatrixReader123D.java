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
	private static int ALPHA = 0, BETA = 1, OTHER = 2;

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
	
	/**
	 * reads a file containing the weights for 123D-Alignment in the following specific manner
	 * 
	 * #Alpha	Beta	Other
	 * 0.0d		0.0d	0.0d	#go
	 * 0.0d		0.0d	0.0d	#ge
	 * 0.0d		0.0d	0.0d	#seq
	 * 0.0d		0.0d	0.0d
	 * 
	 * @param filename
	 * @return the secondary structure preference as a double array
	 */
	public double[][] readWeights(String filename) {
		int[] ind = new int[3];
		
		double[] seq = new double[3];
		double[] go = new double[3];
		double[] ge = new double[3];
		double[] ssp = new double[3];
		double[] lccp = new double[3];
		double[] gccp = new double[3];
		
		double[][] multiresult = new double[6][3];
		
		try {
			BufferedReader in = new BufferedReader(new FileReader(filename));
			String line;
			while ((line = in.readLine()) != null) {
				if(line.startsWith("#")) {
					ind = calcIndices(line.split("\\s+"));
				}
				else {
					String[] cont = line.split("\\s+");
					if(cont[cont.length].equals("go")) {
						go[ALPHA] = Double.parseDouble(cont[ind[ALPHA]]);
						go[BETA] = Double.parseDouble(cont[ind[BETA]]);
						go[OTHER] = Double.parseDouble(cont[ind[OTHER]]);
					}
					else if(cont[cont.length].equals("ge"))
						ge[ALPHA] = Double.parseDouble(cont[ind[ALPHA]]);
					else if(cont[cont.length].equals("seq"))
						seq[ALPHA] = Double.parseDouble(cont[ind[ALPHA]]);
						
				}
			}
			in.close();
		} catch (IOException e) {
			System.out.println("No Input (secondary structure preference)!");
		}
		multiresult[0] = seq;
		multiresult[1] = go;
		multiresult[2] = ge;
		multiresult[3] = ssp;
		multiresult[4] = lccp;
		multiresult[5] = gccp;
		return multiresult;
	}


	private int[] calcIndices(String[] split) {
		int[] result = new int[3];
		for(int i = 1; i < split.length; i++) {
			if(split[i].equals("l"))
				result[i-1] = OTHER;
			if(split[i].equals("h"))
				result[i-1] = ALPHA;
			if(split[i].equals("s"))
				result[i-1] = BETA;
		}
		return result;
	}

}
