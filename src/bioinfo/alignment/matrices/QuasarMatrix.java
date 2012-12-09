package bioinfo.alignment.matrices;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;

public class QuasarMatrix {
	
	public final static double[][] DAYHOFF_MATRIX =
		{
			{ 1.80, -1.60,   .20,   .30, -2.00, - .40,   .30,  1.30, -1.40, - .50, -1.90, -1.20, -1.20, -3.50,  1.10,  1.10,  1.20, -5.60, -3.50,   .20},
			{-1.60,  6.10, -0.00, -1.30, -3.60,  1.20, -1.10, -2.60,  1.50, -2.00, -3.00,  3.40, - .50, -4.50, - .20, - .30, - .90,  2.30, -4.20, -2.50},
			{  .20, -0.00,  2.00,  2.10, -3.60,   .80,  1.40,   .40,  1.60, -1.80, -2.90,  1.00, -1.80, -3.50, - .50,   .70,   .40, -3.90, -2.10, -1.80},
			{  .30, -1.30,  2.10,  3.90, -5.10,  1.60,  3.40,   .60,   .70, -2.40, -4.00,   .10, -2.60, -5.60, -1.00,   .30, - .10, -6.60, -4.30, -2.20},
			{-2.00, -3.60, -3.60, -5.10, 12.00, -5.30, -5.30, -3.30, -3.40, -2.30, -6.00, -5.40, -5.20, -4.30, -2.70,  0.00, -2.20, -7.50,   .40, -1.90},
			{- .40,  1.20,   .80,  1.60, -5.30,  4.10,  2.50, -1.20,  2.90, -2.00, -1.80,   .70, -1.00, -4.70,   .20, - .50, - .80, -4.60, -4.00, -1.90},
			{  .30, -1.10,  1.40,  3.40, -5.30,  2.50,  3.90,   .20,   .60, -2.00, -3.30, - .10, -2.20, -5.40, - .60,  0.00, - .40, -6.80, -4.30, -1.80},
			{ 1.30, -2.60,   .40,   .60, -3.30, -1.20,   .20,  4.80, -2.10, -2.60, -4.00, -1.70, -2.80, -4.80, - .50,  1.10,  0.00, -6.80, -5.30, -1.40},
			{-1.40,  1.50,  1.60,   .70, -3.40,  2.90,   .60, -2.10,  6.60, -2.50, -2.10, - .10, -2.20, -1.80, - .30, - .80, -1.30, -2.50, - .10, -2.30},
			{- .50, -2.00, -1.80, -2.40, -2.30, -2.00, -2.00, -2.60, -2.50,  4.60,  2.40, -1.90,  2.20,  1.00, -2.00, -1.40,   .10, -5.00, -1.00,  3.70},
			{-1.90, -3.00, -2.90, -4.00, -6.00, -1.80, -3.30, -4.00, -2.10,  2.40,  6.00, -2.90,  3.70,  1.80, -2.60, -2.80, -1.70, -1.70, - .90,  1.80},
			{-1.20,  3.40,  1.00,   .10, -5.40,   .70, - .10, -1.70, - .10, -1.90, -2.90,  4.70,   .40, -5.30, -1.20, - .20,  0.00, -3.30, -4.50, -2.50},
			{-1.20, - .50, -1.80, -2.60, -5.20, -1.00, -2.20, -2.80, -2.20,  2.20,  3.70,   .40,  6.60,   .20, -2.10, -1.60, - .60, -4.10, -2.50,  1.80},
			{-3.50, -4.50, -3.50, -5.60, -4.30, -4.70, -5.40, -4.80, -1.80,  1.00,  1.80, -5.30,   .20,  9.10, -4.60, -3.20, -3.10,   .50,  7.00, -1.20},
			{ 1.10, - .20, - .50, -1.00, -2.70,   .20, - .60, - .50, - .30, -2.00, -2.60, -1.20, -2.10, -4.60,  5.90,   .90,   .30, -5.50, -5.00, -1.20},
			{ 1.10, - .30,   .70,   .30, - .00, - .50, - .00,  1.10, - .80, -1.40, -2.80, - .20, -1.60, -3.20,   .90,  1.60,  1.30, -2.30, -2.80, -1.00},
			{ 1.20, - .90,   .40, - .10, -2.20, - .80, - .40, - .00, -1.30,   .10, -1.70, -0.00, - .60, -3.10,   .30,  1.30,  2.60, -5.00, -2.80,   .30},
			{-5.60,  2.30, -3.90, -6.60, -7.50, -4.60, -6.80, -6.80, -2.50, -5.00, -1.70, -3.30, -4.10,   .50, -5.50, -2.30, -5.00, 17.30,  0.00, -6.10},
			{-3.50, -4.20, -2.10, -4.30,   .40, -4.00, -4.30, -5.30, - .10, -1.00, - .90, -4.50, -2.50,  7.00, -5.00, -2.80, -2.80,  0.00, 10.20, -2.50},
			{  .20, -2.50, -1.80, -2.20, -1.90, -1.90, -1.80, -1.40, -2.30,  3.70,  1.80, -2.50,  1.80, -1.20, -1.20, -1.00,   .30, -6.10, -2.50,  4.30}
		};
	
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
