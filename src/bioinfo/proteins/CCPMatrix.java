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

	/**
	 * This is the data in the matrix: first dimension is amino acid [0..25]
	 * second dimension is secondary structure [0..2] third dimension is number
	 * of local contacts [0..13] fourth dimension is number of long-range
	 * contacts [0..13]
	 */
	private double[][][][] matrix;

	public CCPMatrix(double[][][][] matrix) {
		this.matrix = matrix;
	}

	/**
	 * constructs a new CCPMatrix with the Matrix given in the folder. FileNames
	 * must be ccpa ccpla ccpb ccplb ccpo ccplo and ssp (ccpa = global alpha,
	 * ccpb = global beta, ccpo = global coil ccpla = local alpha, ccplb = local
	 * beta, ccplo = local coil ssp = secondaryStructurePotential)
	 * 
	 * @param foldername
	 */
	public CCPMatrix(String foldername) {
		this(readFromFiles(foldername));
	}

	/**
	 * 
	 * @param aa
	 *            AMinoAcidType
	 * @param ss
	 * @param local
	 * @param global
	 * @return the value from the matrix
	 */
	public double getValue(AminoAcidName aa, SecStructThree ss, int local,
			int global) {
		if (ss == SecStructThree.H) { // alpha Helix
			return matrix[aa.getNumber()][0][local][global];
		} else if (ss == SecStructThree.E) { // beta shEEt
			return matrix[aa.getNumber()][1][local][global];
		} else { // Coil
			return matrix[aa.getNumber()][2][local][global];
		}
	}

	/**
	 * 
	 * @param foldername
	 * @return
	 */
	private static double[][][][] readFromFiles(String foldername) {

		double[][][][] result = new double[26][3][14][14];

		if (!foldername.endsWith("/")) {
			foldername = foldername + "/";
		}

		BufferedReader br = null;

		// read ccpa file
		try {
			br = new BufferedReader(new FileReader(foldername + "ccpa"));
			String line = null;
			int aa = 0;
			while ((line = br.readLine()) != null) {
				if (line.trim().startsWith("0")) {
					continue; // skip
				} else {
					aa = AminoAcidName.getNumberFromTLC(line.substring(0, 3));
					for (int i = 0; i <= 13; i++) {
						// TODO
						// matrix[aa][0][][];
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
		return null;
	}

	/**
	 * 
	 * @param foldername
	 * @return
	 */
	private static void writeToFiles(double[][][][] matrix, String foldername) {

		if (!foldername.endsWith("/")) {
			foldername = foldername + "/";
		}

		BufferedWriter bw = null;

		// read ccpa file
		try {
			bw = new BufferedWriter(new FileWriter(foldername + "ccpa"));

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
