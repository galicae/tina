/******************************************************************************
 * bioinfo.energy.potential.hydrophobicity.HydrophobicityMatrix.java          *
 *                                                                            *
 * This class represents an HydrophobicityMatrix.                             *
 *                                                                            *
 * The idea is: w[i][j] = r[i][j]/(p[i]*q[j])                                 *
 * where r[i][j] is the relative frequency of amino acid [i] at dob[j],       *
 *       p[i] is the relative frequency of amino acid [i] and                 *
 *       q[j] is the relative frequency of dob[j]                             *
 *                                                                            *
 * @see Protein Threading by Recursive Dynamic Programming. JMB 290, 757-779  *
 *                                                                            *
 * This file is best read at line width 80 and tab width 4.                   *
 *                                                                   huberste *
 ******************************************************************************/
package bioinfo.energy.potential.hydrophobicity;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.Locale;

/**
 * @author huberste
 * @lastchange 2013-02-17
 */
public class HydrophobicityMatrix {

	// TODO fill matrix with values!
	private final static double[][] STANDARD_MATRIX = {};
	
	private final double[][] matrix;

	/**
	 * constructs a HydrophobicityMatrix
	 * @param matrix
	 */
	public HydrophobicityMatrix(double[][] matrix) {
		this.matrix = matrix;
	}
	
	/**
	 * constructs a HydrophobicityMatrix with standard parameters
	 */
	public HydrophobicityMatrix() {
		this(STANDARD_MATRIX);
	}
	
	/**
	 * constructs a HydrophobicityMatrix with the same parameters as the given
	 * HydrophobicityMatrix
	 * @param arg
	 */
	public HydrophobicityMatrix(HydrophobicityMatrix arg) {
		this(arg.matrix);
	}
	
	/**
	 * constructs a new HydrophobicityMatrix with the Matrix given in the file.
	 * @param filename
	 */
	public HydrophobicityMatrix(String filename) {
		this(readFile(filename));
	}
	
	/**
	 * 
	 * @param aa
	 * @param bucket
	 * @return the value from the matrix
	 */
	public double getValue(int aa, int bucket) {
		return matrix[aa][bucket];
	}
	
	public int getBuckets() {
		return matrix[0].length;
	}
	
	/**
	 * writes a HydrophobicityMatrix file
	 * @param w the HydrophobicityMatrix to be written
	 */
	public static void writeFile(double[][] w, String filepath) {
		
		// initialize important stuff	
		Locale.setDefault(Locale.US);
		DecimalFormat df = new DecimalFormat("0.000000");
		
		BufferedWriter bw = null;
		try {
			bw = new BufferedWriter(new FileWriter(filepath));
			bw.write("# Hydrophobicity matrix\n");
			bw.write("# degree of burial x amino acid\n");
			bw.write("# dob buckets: "+ w[0].length + "\n");
			for (int i = 0; i < w.length; i++) {
				char x = (char)(65+i);
				bw.write(x);
				for (int j = 0; j < w[i].length; j++) {
					if (Double.isNaN(w[i][j])) {
						bw.write("\tNAN");
					} else if(Double.isInfinite(w[i][j])) {
						bw.write("\tNAN");
					} else { 
						bw.write("\t" + df.format(w[i][j]));
					}
				}
				bw.write("\n");
			}
		} catch (IOException e) {
			System.err.println("Error 187: problems reading an outfile:");
			e.printStackTrace();
		} finally {
			try {
				if (bw != null) {
					bw.close();
					bw = null; // give to GarbageCollector
				}
			} catch (IOException e) {
				System.err.println("Error 193: problems closing an outfile:");
				e.printStackTrace();
			}
		}
	}
	
	/**
	 * reads a HydrophobicityMatrix file
	 * @param filepath
	 * @return the double[][] matrix represented by the file
	 */
	public static double[][] readFile(String filepath) {
		BufferedReader br = null;
		double[][] result = null;
		try {
			br = new BufferedReader(new FileReader(filepath));
			String line = null;
			int lines = 0;
			while ((line = br.readLine()) != null) {
				if (line.startsWith("#")) { // comment
					continue;
				}
				if (!line.isEmpty()) // empty lines
					lines++;
			}
			br.close();
			result = new double[lines][];
			lines = 0;
			br = new BufferedReader(new FileReader(filepath));
			while ((line = br.readLine()) != null) {
				if (line.startsWith("#")) { // comment
					continue;
				}
				String[] temp = line.split("\t");
				result[lines] = new double[temp.length-1];
				for (int i = 0; i < temp.length-1; i++) {
					if (temp[i+1].startsWith("NAN")) {
						result[lines][i] = Double.NEGATIVE_INFINITY;
					} else {
						result[lines][i] = Double.parseDouble(temp[i+1]);
					}
				}
				lines++;
			}
			
		} catch (IOException e) {
			System.err.println("Error 221: problems reading the inputfile:");
			e.printStackTrace();
		} finally {
			try {
				if(br != null) {
					br.close();
					br = null;
				}
			} catch (IOException e) {
				System.err.println("Error 230: problems closing the inputfile:");
				e.printStackTrace();
			}
		}
		return result;
	}
	
}

/******************************************************************************
 * "A question that sometimes drives me hazy:                                 *
 *  Am I or are the others crazy?"                                            *
 *     - Albert Einstein (1879 - 1955)                                        *
 ******************************************************************************/
