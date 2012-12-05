package bioinfo.alignment.matrices;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;

/**
 * Reader class for all matrices that 123D needs apart from subsitution matrix
 * @author gobi_4
 * @date December 5, 2012
 *
 */
public class MatrixReader123D {
	public double[][] readSecStructPref(String filename) {
		double [][] matrix = new double[3][21];
		try {
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

}
