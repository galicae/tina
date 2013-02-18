/******************************************************************************
 * bioinfo.proteins.CCPMatrix.java                                            *
 *                                                                            *
 * This file contains the class CCPMatrix which is a representation for       *
 * (Secondary Structure dependant) Contact Capacity Potential Matrices.       *
 *                                                                            *
 * This file is best read at line width 80 and tab width 4.                   *
 *                                                                   huberste *
 ******************************************************************************/
package bioinfo.proteins;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;

/**
 * 
 * @author huberste
 * @lastchange 2013-02-18
 */
public class CCPMatrix {

	private final static String FILE_HEADER = "        0     1     2     3     4     5     6     7     8     9    10    11    12    13\n";

	/**
	 * This is the data in the matrix: first dimension is amino acid [0..25]
	 * second dimension is secondary structure [0..2] third dimension is
	 * local/global third dimension is number of contacts [0..13]
	 */
	private double[][][][] matrix;

	public CCPMatrix(double[][][][] matrix) {
		this.matrix = matrix;
	}

	/**
	 * constructs a new CCPMatrix with the Matrix given in the folder. FileNames
	 * must be $FILENAMEa, $FILENAMEb, $FILENAMEo a = alpha b = beta, o = coil
	 * 
	 * @param foldername
	 */
	public CCPMatrix(String filename) {
		this(readFromFiles(filename));
	}

	/**
	 * 
	 * @param aa
	 *            AMinoAcidType
	 * @param ss
	 *            secondary structure, alpha=0 beta=1 coil=2
	 * @param dist
	 *            local = 0, long-range = 1
	 * @param contacts
	 *            number of contacts
	 * @return the value from the matrix
	 */
	public double getValue(AminoAcidName aa, SecStructThree ss, int dist,
			int contacts) {
		if (ss == SecStructThree.H) { // alpha Helix
			return matrix[aa.getNumber()][0][dist][contacts];
		} else if (ss == SecStructThree.E) { // beta shEEt
			return matrix[aa.getNumber()][1][dist][contacts];
		} else { // Coil
			return matrix[aa.getNumber()][2][dist][contacts];
		}
	}
	
	/**
	 * 
	 * @param aa
	 * @param ss
	 * @param dist
	 * @param contacts
	 * @return
	 */
	public double getValue(char aa, SecStructThree ss, int dist,
			int contacts) {
		if (ss == SecStructThree.H) { // alpha Helix
			return matrix[aa-65][0][dist][contacts];
		} else if (ss == SecStructThree.E) { // beta shEEt
			return matrix[aa-65][1][dist][contacts];
		} else { // Coil
			return matrix[aa-65][2][dist][contacts];
		}
	}

	/**
	 * 
	 * @param filename Filename of the CCP files, not includin a,la, b,lb, o,lo
	 * @return
	 */
	private static double[][][][] readFromFiles(String filename) {

		double[][][][] result = new double[26][3][2][14];

		double[][] resultccpa = readCCPFile(filename + "a", 0);
		double[][] resultccpb = readCCPFile(filename + "b", 1);
		double[][] resultccpo = readCCPFile(filename + "o", 2);
		double[][] resultccpla = readCCPFile(filename + "la", 0);
		double[][] resultccplb = readCCPFile(filename + "lb", 1);
		double[][] resultccplo = readCCPFile(filename + "lo", 2);

		for (int aa = 0; aa < 26; aa++) {
			result[aa][0][0] = resultccpla[aa];
			result[aa][1][0] = resultccplb[aa];
			result[aa][2][0] = resultccplo[aa];
			result[aa][0][1] = resultccpa[aa];
			result[aa][1][1] = resultccpb[aa];
			result[aa][2][1] = resultccpo[aa];
			
		}

		return result;
	}

	private static double[][] readCCPFile(String filename, int ss) {

		double[][] result = new double[26][14];

		BufferedReader br = null;
		try {
			br = new BufferedReader(new FileReader(filename));
			String line = null;
			int aa = 0;
			while ((line = br.readLine()) != null) {
				if (line.trim().startsWith("0")) {
					continue; // skip
				} else {
					aa = AminoAcidName.getNumberFromTLC(line.substring(0, 3));
					for (int i = 0; i <= 13; i++) {
						result[aa][i] = Double.parseDouble(line.substring(
								3 + (i * 6), 9 + (i * 6)).trim());
					}
				}

			}
		} catch (IOException e) {
			System.err.println("Error occured while reading a file:"
					+ e.getLocalizedMessage());
			e.printStackTrace();
		} finally {
			try {
				if (br != null) {
					br.close(); // RessourceLeak
					br = null; // set null -> GC
				}
			} catch (IOException e) {
				System.err.println("Error occured while closing a file:"
						+ e.getLocalizedMessage());
				e.printStackTrace();
			}
		}
		return result;
	}

	/**
	 * 
	 * @param foldername
	 * @return
	 */
	public static void writeToFiles(double[][][][] matrix, String filename) {
		writeCCPFile(filename + "la", 0, 0, matrix);
		writeCCPFile(filename + "lb", 1, 0, matrix);
		writeCCPFile(filename + "lo", 2, 0, matrix);
		writeCCPFile(filename +  "a", 0, 1, matrix);
		writeCCPFile(filename +  "b", 1, 1, matrix);
		writeCCPFile(filename +  "o", 2, 1, matrix);
	}

	private static void writeCCPFile(String filename, int ss, int dist,
			double[][][][] matrix) {
		BufferedWriter bw = null;
		try {
			bw = new BufferedWriter(new FileWriter(filename + "b"));
			bw.write(FILE_HEADER);
			for (int aa = 0; aa < 26; aa++) {
				bw.write(AminoAcidName.getTLCFromNumber(aa));
				for (int contacts = 0; contacts <= 13; contacts++) {
					String write = Double.toString(matrix[aa][ss][dist][contacts]);
					while (write.length() < 6)
						// needs to be length 6
						write = " " + write;
					bw.write(write);
				}
				bw.write("\n");
			}

		} catch (IOException e) {
			System.err.println("Error occured while reading a file:"
					+ e.getLocalizedMessage());
			e.printStackTrace();
		} finally {
			try {
				if (bw != null) {
					bw.close(); // RessourceLeak
					bw = null; // set null -> GC
				}
			} catch (IOException e) {
				System.err.println("Error occured while closing a file:"
						+ e.getLocalizedMessage());
				e.printStackTrace();
			}
		}
	}

	/**
	 * for testing
	 * 
	 * @param args
	 */
	public static void main(String[] args) {

	}

}

/******************************************************************************
 * "A question that sometimes drives me hazy: * Am I or are the others crazy?" *
 * - Albert Einstein (1879 - 1955) *
 ******************************************************************************/
