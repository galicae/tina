/******************************************************************************
 * util.MatrixTools.java                                                      *
 *                                                                            *
 * This class provides some tools for working with matrices.                  *
 *                                                                            *
 * This file is best read at line width 80 and tab width 4.                   *
 *                                                                   huberste *
 ******************************************************************************/
package util;

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
	
}
