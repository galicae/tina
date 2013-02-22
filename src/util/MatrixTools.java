/******************************************************************************
 * util.MatrixTools.java                                                      *
 *                                                                            *
 * This class provides some tools for working with matrices.                  *
 *                                                                            *
 * This file is best read at line width 80 and tab width 4.                   *
 *                                                                   huberste *
 ******************************************************************************/
package util;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Locale;

/**
 * @author huberste
 *
 */
public class MatrixTools {

	public static String toJavaString(double[][] matrix) {
		String result = "{";
		for (int x = 0; x < matrix.length-1; x++) {
			result += toJavaString(matrix[x]) + ", ";
		}
		
		result +="}";
		return result;
	}
	
	/**
	 * returns a vector as a Java Matrix String
	 * @param matrix
	 * @return
	 */
	public static String toJavaString(double[] matrix) {
		String result = "{";
		for (int x = 0; x < matrix.length-1; x++) {
			result += matrix[x] + ", ";
		}
		result += matrix[matrix.length -1];
		result +="}";
		return result;
	}

	/**
	 * TEST
	 * @param objectarray
	 * @return
	 */
	public static String toJavaStringx(Object[] objectarray) {
		String result = "{";
		for (int x = 0; x < objectarray.length-1; x++) {
			result += objectarray[x] + ", ";
		}
		result += objectarray[objectarray.length -1];
		result +="}";
		return result;
	}
	
	public static void writeToFile(double[][] matrix, String fileName, int precision) {
		BufferedWriter bw = null;
		try {
			bw = new BufferedWriter(new FileWriter(fileName));
			int y = 0;
			while (y < matrix[0].length) {
				int x = 0;
				while(x < matrix.length) {
					bw.write(String.format(Locale.US,"%."+precision+"f", matrix[x][y])+"\t");
					x++;
				}
				bw.write("\n");
				y++;
			}
		} catch (IOException e) {
			System.err.println("Error writing file: " + e.getLocalizedMessage());
			e.printStackTrace();
		} finally {
			try {
				if (bw != null) {
					bw.close();
					bw = null;	// GC
				}
			} catch (IOException e) {
				System.err.println("Error closing file: " + e.getLocalizedMessage());
				e.printStackTrace();
			}
		}
	}
	
}
