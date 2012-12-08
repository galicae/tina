package bioinfo.alignment.matrices;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;

public class QuasarMatrix {
	
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
