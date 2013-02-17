/******************************************************************************
 * bioinfo.energy.potential.contactcapacity.ContactCapacityMatrix.java        *
 *                                                                            *
 * This class represents an ContactCapacityMatrix.                            *
 *                                                                            *
 * @see Fast Protein Fold Recognition via Sequence to Structure Alignment and *
 * Contact Capacity Potentials. Pacific Symposium on Biocomputing 96, pp53-72 *
 *                                                                            *
 * This file is best read at line width 80 and tab width 4.                   *
 *                                                                   huberste *
 ******************************************************************************/
package bioinfo.energy.potential.contactcapacity;

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
public class ContactCapacityMatrix {

	// TODO extend to secondary structure dependent ContactCapacityPotential or
	//	ConditionalContactCapacityPotential
	
	/**
	 * maximum residue distance for a contact to be considered local
	 */
	public final static int MAX_LOCAL_DISTANCE = 4;
	
	/**
	 * Maximum distance between two (virtual) C-beta Atoms that defines a
	 * contact (in Angstrom)
	 */
	public final static double MAX_CONTACT_DISTANCE = 7.0;
	
	/**
	 * this is the number of buckets for local contacts
	 */
	public final static int MAX_LOCAL_CONTACTS= 10;
	
	/**
	 * this is the number of buckets fot long-range contacts
	 */
	public final static int MAX_LONGRANGE_CONTACTS=12;
	
	// TODO fill matrix with values!
	private final static double[][][] STANDARD_MATRIX = {};
	
	/**
	 * matrix[AA][local][longrange]
	 */
	private final double[][][] matrix;

	/**
	 * constructs a ContactCapacityMatrix
	 * @param matrix
	 */
	public ContactCapacityMatrix(double[][][] matrix) {
		this.matrix = matrix;
	}
	
	/**
	 * constructs a ContactCapacityMatrix with standard parameters
	 */
	public ContactCapacityMatrix() {
		this(STANDARD_MATRIX);
	}
	
	/**
	 * constructs a ContactCapacityMatrix with the same parameters as the given
	 * HydrophobicityMatrix
	 * @param arg
	 */
	public ContactCapacityMatrix(ContactCapacityMatrix arg) {
		this(arg.matrix);
	}
	
	/**
	 * constructs a new ContactCapacityMatrix with the Matrix given in the file.
	 * @param filename
	 */
	public ContactCapacityMatrix(String filename) {
		this(readCCPMatrixFromFile(filename));
	}
	
	/**
	 * 
	 * @param aa
	 * @param localContacts
	 * @param longRangeContacts
	 * @return the value from the matrix
	 */
	public double getValue(int aa, int localContacts, int longRangeContacts) {
		return matrix[aa][localContacts][longRangeContacts];
	}
	
	/**
	 * writes a ContactCapacityMatrix file
	 * @param w the ContactCapacityMatrix to be written
	 */
	public static void writeCCPMatrixToFile(double[][][] w, String filepath) {
		// TODO
		
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
					bw.write("\t" + df.format(w[i][j]));
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
	 * reads a ContactCapacityMatrix file
	 * @param filepath
	 * @return the double[][][] matrix represented by the file
	 */
	public static double[][][] readCCPMatrixFromFile(String filepath) {
		// TODO
		BufferedReader br = null;
		double[][][] result = null;
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
			result = new double[lines][][];
			lines = 0;
			br = new BufferedReader(new FileReader(filepath));
			while ((line = br.readLine()) != null) {
				if (line.startsWith("#")) { // comment
					continue;
				}
				String[] temp = line.split("\t");
//				result[lines] = new double[temp.length-1];
				for (int i = 0; i < temp.length-1; i++) {
//					result[lines][i] = Double.parseDouble(temp[i+1]);
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
