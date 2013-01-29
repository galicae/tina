package bioinfo.alignment.matrices;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;

public class QuasarMatrix {
	
	/**
	 * hard-coded dayhoff matrix. Hark at this. Proves really useful for testing since it is a default matrix.
	 */
	public final static double[][] DAYHOFF_MATRIX =
		{
			{1.8, 0.0, -2.0, 0.3, 0.3, -3.5, 1.3, -1.4, -0.5, 0.0, -1.2, -1.9, -1.2, 0.2, 0.0, 1.1, -0.4, -1.6, 1.1, 1.2, 0.0, 0.2, -5.6, 0.0, -3.5, 0.0},
			{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
			{-2.0, 0.0, 12.0, -5.1, -5.3, -4.3, -3.3, -3.4, -2.3, 0.0, -5.4, -6.0, -5.2, -3.6, 0.0, -2.7, -5.3, -3.6, -0.0, -2.2, 0.0, -1.9, -7.5, 0.0, 0.4, 0.0},
			{0.3, 0.0, -5.1, 3.9, 3.4, -5.6, 0.6, 0.7, -2.4, 0.0, 0.1, -4.0, -2.6, 2.1, 0.0, -1.0, 1.6, -1.3, 0.3, -0.1, 0.0, -2.2, -6.6, 0.0, -4.3, 0.0},
			{0.3, 0.0, -5.3, 3.4, 3.9, -5.4, 0.2, 0.6, -2.0, 0.0, -0.1, -3.3, -2.2, 1.4, 0.0, -0.6, 2.5, -1.1, -0.0, -0.4, 0.0, -1.8, -6.8, 0.0, -4.3, 0.0},
			{-3.5, 0.0, -4.3, -5.6, -5.4, 9.1, -4.8, -1.8, 1.0, 0.0, -5.3, 1.8, 0.2, -3.5, 0.0, -4.6, -4.7, -4.5, -3.2, -3.1, 0.0, -1.2, 0.5, 0.0, 7.0, 0.0},
			{1.3, 0.0, -3.3, 0.6, 0.2, -4.8, 4.8, -2.1, -2.6, 0.0, -1.7, -4.0, -2.8, 0.4, 0.0, -0.5, -1.2, -2.6, 1.1, -0.0, 0.0, -1.4, -6.8, 0.0, -5.3, 0.0},
			{-1.4, 0.0, -3.4, 0.7, 0.6, -1.8, -2.1, 6.6, -2.5, 0.0, -0.1, -2.1, -2.2, 1.6, 0.0, -0.3, 2.9, 1.5, -0.8, -1.3, 0.0, -2.3, -2.5, 0.0, -0.1, 0.0},
			{-0.5, 0.0, -2.3, -2.4, -2.0, 1.0, -2.6, -2.5, 4.6, 0.0, -1.9, 2.4, 2.2, -1.8, 0.0, -2.0, -2.0, -2.0, -1.4, 0.1, 0.0, 3.7, -5.0, 0.0, -1.0, 0.0},
			{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
			{-1.2, 0.0, -5.4, 0.1, -0.1, -5.3, -1.7, -0.1, -1.9, 0.0, 4.7, -2.9, 0.4, 1.0, 0.0, -1.2, 0.7, 3.4, -0.2, -0.0, 0.0, -2.5, -3.3, 0.0, -4.5, 0.0},
			{-1.9, 0.0, -6.0, -4.0, -3.3, 1.8, -4.0, -2.1, 2.4, 0.0, -2.9, 6.0, 3.7, -2.9, 0.0, -2.6, -1.8, -3.0, -2.8, -1.7, 0.0, 1.8, -1.7, 0.0, -0.9, 0.0},
			{-1.2, 0.0, -5.2, -2.6, -2.2, 0.2, -2.8, -2.2, 2.2, 0.0, 0.4, 3.7, 6.6, -1.8, 0.0, -2.1, -1.0, -0.5, -1.6, -0.6, 0.0, 1.8, -4.1, 0.0, -2.5, 0.0},
			{0.2, 0.0, -3.6, 2.1, 1.4, -3.5, 0.4, 1.6, -1.8, 0.0, 1.0, -2.9, -1.8, 2.0, 0.0, -0.5, 0.8, -0.0, 0.7, 0.4, 0.0, -1.8, -3.9, 0.0, -2.1, 0.0},
			{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
			{1.1, 0.0, -2.7, -1.0, -0.6, -4.6, -0.5, -0.3, -2.0, 0.0, -1.2, -2.6, -2.1, -0.5, 0.0, 5.9, 0.2, -0.2, 0.9, 0.3, 0.0, -1.2, -5.5, 0.0, -5.0, 0.0},
			{-0.4, 0.0, -5.3, 1.6, 2.5, -4.7, -1.2, 2.9, -2.0, 0.0, 0.7, -1.8, -1.0, 0.8, 0.0, 0.2, 4.1, 1.2, -0.5, -0.8, 0.0, -1.9, -4.6, 0.0, -4.0, 0.0},
			{-1.6, 0.0, -3.6, -1.3, -1.1, -4.5, -2.6, 1.5, -2.0, 0.0, 3.4, -3.0, -0.5, -0.0, 0.0, -0.2, 1.2, 6.1, -0.3, -0.9, 0.0, -2.5, 2.3, 0.0, -4.2, 0.0},
			{1.1, 0.0, -0.0, 0.3, -0.0, -3.2, 1.1, -0.8, -1.4, 0.0, -0.2, -2.8, -1.6, 0.7, 0.0, 0.9, -0.5, -0.3, 1.6, 1.3, 0.0, -1.0, -2.3, 0.0, -2.8, 0.0},
			{1.2, 0.0, -2.2, -0.1, -0.4, -3.1, -0.0, -1.3, 0.1, 0.0, -0.0, -1.7, -0.6, 0.4, 0.0, 0.3, -0.8, -0.9, 1.3, 2.6, 0.0, 0.3, -5.0, 0.0, -2.8, 0.0},
			{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
			{0.2, 0.0, -1.9, -2.2, -1.8, -1.2, -1.4, -2.3, 3.7, 0.0, -2.5, 1.8, 1.8, -1.8, 0.0, -1.2, -1.9, -2.5, -1.0, 0.3, 0.0, 4.3, -6.1, 0.0, -2.5, 0.0},
			{-5.6, 0.0, -7.5, -6.6, -6.8, 0.5, -6.8, -2.5, -5.0, 0.0, -3.3, -1.7, -4.1, -3.9, 0.0, -5.5, -4.6, 2.3, -2.3, -5.0, 0.0, -6.1, 17.3, 0.0, 0.0, 0.0},
			{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
			{-3.5, 0.0, 0.4, -4.3, -4.3, 7.0, -5.3, -0.1, -1.0, 0.0, -4.5, -0.9, -2.5, -2.1, 0.0, -5.0, -4.0, -4.2, -2.8, -2.8, 0.0, -2.5, 0.0, 0.0, 10.2, 0.0},
			{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}
		};
	
	/**
	 * parses a substitution matrix. Expects to see a file formatted as a quasar matrix.
	 * @param filename the matrix input file
	 * @return the matrix as a double[][] array
	 */
	public static double[][] parseMatrix(String filename){
		double [][] matrix = new double[26][26];
		try{
			BufferedReader in = new BufferedReader(new FileReader(filename));
			String line;
			int numrows = 0;
			char[] aminos = null; 
			while((line = in.readLine()) != null){
				if (line.matches("^(NUMROW)(.)*")){
					numrows = Integer.parseInt(line.substring(7,9));
				}
				else if (line.matches("^(ROWINDEX)(.)*")){
					aminos = line.substring(9,numrows+9).toCharArray();
				}
				else if(line.matches("^(MATRIX)(.)*")){
					String temp[];
					for (int i = 0; i < aminos.length; i++) {
						temp = line.split("\\s+");
						for (int j = 1; j < temp.length; j++) {
							matrix[aminos[i]-65][aminos[j-1]-65] = Double.parseDouble(temp[j]);
							matrix[aminos[j-1]-65][aminos[i]-65] = Double.parseDouble(temp[j]);
						}
						line = in.readLine();
					}
				}
			}
			in.close();
		} catch(IOException e){
			System.out.println("No Input (quasarmatrix)!");
		}
		return matrix;
	}
	
	/**
	 * 
	 * @param matrixString
	 * @param arg
	 * @return
	 */
	public static double[][] parseMatrix(String matrixString, boolean arg){
		double [][] matrix = new double[26][26];
		String[] lines = matrixString.split("\n");
		try{
			String line;
			int currline = 0;
			int numrows = 0;
			char[] aminos = null; 
			while(currline < lines.length){
				line=lines[currline];
				if (line.matches("^(NUMROW)(.)*")){
					numrows = Integer.parseInt(line.substring(7,9));
				}
				else if (line.matches("^(ROWINDEX)(.)*")){
					aminos = line.substring(9,numrows+9).toCharArray();
				}
				else if(line.matches("^(MATRIX)(.)*")){
					String temp[];
					for (int i = 0; i < aminos.length; i++) {
						temp = line.split("\\s+");
						for (int j = 1; j < temp.length; j++) {
							matrix[aminos[i]-65][aminos[j-1]-65] = Double.parseDouble(temp[j]);
							matrix[aminos[j-1]-65][aminos[i]-65] = Double.parseDouble(temp[j]);
						}
						currline++;
						// TODO: maybe this can be made more efficient?
						if (currline < lines.length) {
							line=lines[currline];
						}
					}
				}
				currline++;
			}
		} catch (Exception e){
			System.err.println("Apparently a wrong matrix was given! (QuasarMatrix.parseMatrix(String, boolean))");
			e.printStackTrace();
		}
		return matrix;
	}
	
	public static String getDoubleMatrixAsJavaString(double[][] arg) {
		String matrixAsJavaArray = "{\n\t";
		for (int row = 0; row < arg.length; row ++) {
			matrixAsJavaArray += "{";
			for (int col = 0; col < arg[row].length-1; col++) {
				matrixAsJavaArray += String.valueOf(arg[row][col]);
				if (col < arg[row].length -1) matrixAsJavaArray += ", ";
			}
			matrixAsJavaArray += "}";
			if (row < arg.length -1) matrixAsJavaArray += ",\n\t";
			else matrixAsJavaArray += "\n";
		}
		matrixAsJavaArray += "}";
		return matrixAsJavaArray;
	}
	
//	public static void main(String[] args){
//		int[][] matrix = QuasarMatrix.parseMatrix(args[0]);
//		for (int i = 0; i < matrix[0].length; i++) {
//		for (int j = 0; j < matrix.length; j++) {
//			System.out.print(matrix[j][i] + " | ");
//		}
//			System.out.println();
//		}
//	}
}
