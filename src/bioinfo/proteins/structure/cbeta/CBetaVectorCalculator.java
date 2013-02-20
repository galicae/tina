/******************************************************************************
 * bioinfo.proteins.structure.cbeta.CBetaVectorCalculator.java                *
 *                                                                            *
 * This file is best read at line width 80 and tab width 4.                   *
 *                                                                   huberste *
 ******************************************************************************/
package bioinfo.proteins.structure.cbeta;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.LinkedList;
import java.util.Locale;

import bioinfo.pdb.PDBFile;
import bioinfo.proteins.Atom;
import bioinfo.proteins.AtomType;
import bioinfo.proteins.PDBEntry;
import bioinfo.proteins.PDBFileReader;

/**
 * 
 * @author huberste
 * @lastchange 2013-02-17
 */
public class CBetaVectorCalculator {

	public final static String usage = 
			"usage:\n" +
			"\tjava ContactCounter <pdblist> <pdbpath>\n\n" +
			"where <pdblist> is a list of PDB IDs, <pdbpath> is a\n"+
			"(writable) path to the directory containing the PDB files";
	
	public final static int DIMENSIONS = 3;
	
	public static final double[] ORIGIN = { 0.0, 0.0, 0.0 };
	
	public final static String filedesc =
			"# CBeta Vector data file";
	
	/**
	 * @param args
	 */
	public static void main(String[] args) {

		if (args.length < 2) {
			System.out.println(usage);
			System.exit(1);
		}
		String pdbList = args[0];
		String pdbpath = args[1];
		
		// load file
		LinkedList<String> pdbIDs = new LinkedList<String>();
		BufferedReader br = null;
		String line = null;
		
		try {
			br = new BufferedReader(new FileReader(pdbList));
			while ((line = br.readLine()) != null) {
				if (line.trim().startsWith("#")) { // comment
					continue;
				}
				if (!line.trim().isEmpty())	// e.g. last line
					pdbIDs.add(line.trim());
			}
		} catch (IOException e) {
			System.err.println("Error 60: problems reading the File:");
			e.printStackTrace();
		} finally {
			try {
				if (br != null) {
					br.close();
					br = null; // give br to GarbageCollector
				}
			} catch (IOException e){
				System.err.println("Error 69: problems closing the File:");
				e.printStackTrace();
			}
		}
		
		PDBFileReader pdbreader = new PDBFileReader(pdbpath);
		
		LinkedList<double[]>[] vectors1 = new LinkedList[26]; // Calpha -> Cbeta
		LinkedList<double[]>[] vectors2 = new LinkedList[26]; // O -> Calpha
		for (int aa = 0; aa < 26; aa++) {
			vectors1[aa] = new LinkedList<double[]>();
			vectors2[aa] = new LinkedList<double[]>();
		}
		
		int debug = 0;
		for (String id : pdbIDs) { // for each file
			debug++;
			try {
				// begin debugging
				System.out.println("working on id "+id + " ("+debug+" of " + pdbIDs.size()+")");
				// end debugging
				
				PDBEntry structure = pdbreader.readPDBFromFile(
						PDBFile.getFile(pdbpath, id.substring(0,4)),id.charAt(4));
				
				if (structure == null) { // dirty debugging
					System.err.println("Error 88: Error getting ID "+id);
					continue;
				}
				
				// for each amino acid in that pdb
				for (int pos = 0; pos < structure.length(); pos++) {
					int astype = (structure.getAminoAcid(pos).getName().getOneLetterCode().charAt(0))-65;
					// get AminoAcids
					Atom calpha = structure.getAminoAcid(pos).getAtomByType(AtomType.CA);
					Atom cbeta = structure.getAminoAcid(pos).getAtomByType(AtomType.CB);
					Atom o = structure.getAminoAcid(pos).getAtomByType(AtomType.O);
					// generate vectors
					double[] vector1 = new double[3]; // Calpha -> Cbeta
					double[] vector2 = new double[3]; // O -> Calpha
					if (calpha != null && cbeta != null && o != null) {
						for (int dim = 0; dim < DIMENSIONS; dim++) {
							vector1[dim] = cbeta.getPosition()[dim] -
									calpha.getPosition()[dim];
						}
						for (int dim = 0; dim < DIMENSIONS; dim++) {
							vector2[dim] = calpha.getPosition()[dim] -
									o.getPosition()[dim];
						}
						vectors1[astype].add(vector1);
						vectors2[astype].add(vector2);
					}
				} // for each amino acid in that pdb
			} catch (Exception e) { // dirty debugging
				System.err.println("problem at id "+id + " ("+debug+" of " + pdbIDs.size()+")");
				e.printStackTrace();
			}
		} // for each pdb
		
		// Rotate each vector
		for(int aa = 0; aa < 26; aa++) {
			if (vectors1[aa].size() > 1) {
				double[] vector1null = vectors1[aa].get(0);
//				double[] vector2null = vectors2[aa].get(0);
				for (int num = 1; num < vectors1[aa].size(); num++) {
					double[] vector1 = vectors1[aa].get(num);
					double[] cross = crossProduct(vector1null, vector1);
					double scalar = scalarProduct(vector1null, vector1);
					double phi = Math.acos(scalar / (vectorLength(vector1null)*vectorLength(vector1)));
					double[][] rotationMatrix = calcRotationMatrix(cross, phi);
					System.out.println("vector1null:");
					System.out.println(vectorToString(vector1null));
					System.out.println("vector1:");
					System.out.println(vectorToString(vector1));
					
					
					System.out.println("Rotationsmatrix:");
					System.out.println(matrixToString(rotationMatrix));
					
					double[] vector2 = vectors2[aa].get(num);
					double[][] newvector = matrixMultiplication(vector2, rotationMatrix);
					System.out.println("vector2"+"\n"+vectorToString(vector2));
					System.out.println("newvector"+"\n"+matrixToString(newvector));
				}
			}
		}
		
		
		// print out frequencies
		System.out.println(filedesc);
		char x = 0;
		for (int aa = 0; aa < 26; aa++) {
			
		}
			
		
		
		
		// TODO save vector
		// calculate average vector vor glycin

	}
	
	public static double[][] calcRotationMatrix(double[] cross, double phi) {
		double[][] result = {
			{cross[0] * cross[0] * (1-Math.cos(phi)) +            Math.cos(phi),
			 cross[1] * cross[0] * (1-Math.cos(phi)) + cross[2] * Math.sin(phi),
			 cross[2] * cross[0] * (1-Math.cos(phi)) - cross[1] * Math.sin(phi)},
			{cross[0] * cross[1] * (1-Math.cos(phi)) - cross[2] * Math.sin(phi),
			 cross[1] * cross[1] * (1-Math.cos(phi)) +            Math.cos(phi),
			 cross[2] * cross[1] * (1-Math.cos(phi)) + cross[0] * Math.sin(phi)},
			{cross[0] * cross[2] * (1-Math.cos(phi)) + cross[1] * Math.sin(phi),
			 cross[1] * cross[2] * (1-Math.cos(phi)) - cross[0] * Math.sin(phi),
			 cross[2] * cross[2] * (1-Math.cos(phi)) +            Math.cos(phi)}};
		return result;
	}

	public static double vectorLength(double[] vector) {
		return Math.sqrt(Math.pow(vector[0], 2) +
				Math.pow(vector[1], 2) + Math.pow(vector[2], 2));
	}
	
	public static double scalarProduct(double[] vector1, double[] vector2) {
		double result = vector1[0]*vector2[0] +	vector1[1]*vector2[1] +
				vector1[2]*vector2[2]; 
		return result;
	}
	
	public static double[] unitVector(double[] vector) {
		return linearProduct( 1.0 / vectorLength(vector), vector);
	}
	
	public static double[] linearProduct(double factor, double[] vector) {
		double[] result = new double[vector.length];
		for (int i = 0; i < vector.length; i++) {
			result[i] = factor * vector[i];
		}
		return result;
	}
	
	public static double[] crossProduct(double[] vector1, double[] vector2) {
		double[] result =
			{vector1[1]*vector2[2] - vector1[2]*vector2[1],
			 vector1[2]*vector2[0] - vector1[0]*vector2[2],
			 vector1[0]*vector2[1] - vector1[1]*vector2[0]}; 
		return result;
	}
	
	public static double[][] matrixMultiplication(double[][] a, double[][] b) {
		if (a[0].length != b.length) {
			System.err.println("cannot multiply these two matrices!");
			return null;
		}
		
		double[][] result = new double[a.length][b[0].length];
		
		for (int i = 0; i < a.length; i ++) {
			for(int j = 0; j < b[0].length; j++) {
				for (int k = 0; k < a[0].length; k++) {
					result[i][j] += a[i][k] * b[k][j];
				}
			}
		}
		return result;
	}
	
	public static double[][] matrixMultiplication(double[] x, double[][] b) {
		double[][] a = new double[1][];
		a[0] = x;
		if (a[0].length != b.length) {
			System.err.println("cannot multiply these two matrices!");
			return null;
		}
		
		double[][] result = new double[a.length][b[0].length];
		
		for (int i = 0; i < a.length; i ++) {
			for(int j = 0; j < b[0].length; j++) {
				for (int k = 0; k < a[0].length; k++) {
					result[i][j] += a[i][k] * b[k][j];
				}
			}
		}
		return result;
	}
	
	public static String vectorToString(double[] vector) {
		// initialize important stuff	
		Locale.setDefault(Locale.US);
		DecimalFormat df = new DecimalFormat("0.00");
		
		String result = "";
		for (int i = 0; i < vector.length; i++) {
			result += df.format(vector[i]) + "\n";
		}
		return result;
	}
	
	public static String matrixToString(double[][] matrix) {
		// initialize important stuff	
		Locale.setDefault(Locale.US);
		DecimalFormat df = new DecimalFormat("0.00");
		String result = "";
		for(int y = 0; y < matrix[0].length; y++) {
			for (int x = 0; x < matrix.length; x++) {
				result += df.format(matrix[x][y]) + "\t";
			}
			result += "\n";
		}
		return result;
	}

}

/**
 * on the math:
 * n = (n_1, n_2, n_3)^T is rotation unit vector; \alpha is rotationg angle
 * R_{\hat n}(\alpha) = \begin{pmatrix} 
 * n_1^2 \left(1-\cos\alpha\right) + \cos\alpha  & n_1 n_2 \left(1-\cos\alpha\right) - n_3 \sin\alpha &  n_1 n_3 \left(1-\cos\alpha\right) + n_2 \sin\alpha \\
 * n_2 n_1 \left(1-\cos\alpha\right) + n_3 \sin\alpha  & n_2^2\left(1-\cos\alpha\right) + \cos\alpha &   n_2 n_3 \left(1-\cos\alpha\right) - n_1 \sin\alpha         \\
 * n_3 n_1 \left(1-\cos\alpha\right) - n_2 \sin\alpha &  n_3 n_2 \left(1-\cos\alpha\right) + n_1 \sin\alpha & n_3^2\left(1-\cos\alpha\right) + \cos\alpha
 * \end{pmatrix} 
 */

/******************************************************************************
 * "A question that sometimes drives me hazy:                                 *
 *  Am I or are the others crazy?"                                            *
 *     - Albert Einstein (1879 - 1955)                                        *
 ******************************************************************************/
