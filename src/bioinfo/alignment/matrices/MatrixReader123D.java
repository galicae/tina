package bioinfo.alignment.matrices;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;

import util.Bio;

/**
 * Reader class for all matrices that 123D needs apart from subsitution matrix
 * 
 * @author gobi_4
 * @date December 5, 2012
 * 
 */
public class MatrixReader123D {
	private static int ALPHA = 0, BETA = 1, OTHER = 2;

	/**
	 * reads a file containing the secondary structure preferences for all amino
	 * acids given in 3-letter code, with a caps header (ALPHA-BETA-OTHER) over
	 * each section and an UNK entry
	 * 
	 * @param filename
	 * @return the secondary structure preference as a double array
	 */
	public static double[][] readSecStructPref(String filename) {
		double[][] matrix = new double[3][26];
		try {
			BufferedReader in = new BufferedReader(new FileReader(filename));
			String line;
			String[] content = new String[2];
			String secstruct = null;
			char c = 0;
			while ((line = in.readLine()) != null) {			
					content = line.split("\\s+");
					
					if(content[0].equals("ALPHA") || content[0].equals("BETA") || content[0].equals("OTHER")){
						secstruct = content[0];
						line = in.readLine();
						content = line.split("\\s+");
					}
					c = Bio.codeTranslate(content[0]);
					if(secstruct.equals("ALPHA")){	
						matrix[0][c - 65] = Double.parseDouble(content[1]);
					}
					else if(secstruct.equals("BETA")){
						matrix[1][c - 65] = Double.parseDouble(content[1]);
					}
					else if(secstruct.equals("OTHER")){
						matrix[2][c - 65] = Double.parseDouble(content[1]);
					}
			}
			in.close();
		} catch (IOException e) {
			System.out.println("No Input (secondary structure preference)!");
		}
		return matrix;
	}

	public static double[][] readWeights(String filename) {
		int[] ind = new int[3];

		double[] seq = new double[3];
		double[] go = new double[3];
		double[] ge = new double[3];
		double[] ssp = new double[3];
		double[] lccp = new double[3];
		double[] gccp = new double[3];

		double[][] multiresult = new double[6][3];

		String[] cont = new String[6];
		
		try {
			BufferedReader in = new BufferedReader(new FileReader(filename));
			String line;
			while ((line = in.readLine()) != null) {
				cont = line.split("\\s+");
				if (line.startsWith("#")) {
					ind = calcIndices(line.split("\\s+"));
				} else {
					if (cont[cont.length-1].equals("go")) {
						go[ALPHA] = Double.parseDouble(cont[ind[ALPHA]]);
						go[BETA] = Double.parseDouble(cont[ind[BETA]]);
						go[OTHER] = Double.parseDouble(cont[ind[OTHER]]);
					} else if (cont[cont.length-1].equals("ge")) {
						ge[ALPHA] = Double.parseDouble(cont[ind[ALPHA]]);
						ge[BETA] = Double.parseDouble(cont[ind[BETA]]);
						ge[OTHER] = Double.parseDouble(cont[ind[OTHER]]);
					} else if (cont[cont.length-1].equals("seq")) {
						seq[ALPHA] = Double.parseDouble(cont[ind[ALPHA]]);
						seq[BETA] = Double.parseDouble(cont[ind[BETA]]);
						seq[OTHER] = Double.parseDouble(cont[ind[OTHER]]);
					} else if (cont[cont.length-1].equals("ssp")) {
						ssp[ALPHA] = Double.parseDouble(cont[ind[ALPHA]]);
						ssp[BETA] = Double.parseDouble(cont[ind[BETA]]);
						ssp[OTHER] = Double.parseDouble(cont[ind[OTHER]]);
					} else if (cont[cont.length-1].equals("lccp")) {
						lccp[ALPHA] = Double.parseDouble(cont[ind[ALPHA]]);
						lccp[BETA] = Double.parseDouble(cont[ind[BETA]]);
						lccp[OTHER] = Double.parseDouble(cont[ind[OTHER]]);
					} else if (cont[cont.length-1].equals("gccp")) {
						gccp[ALPHA] = Double.parseDouble(cont[ind[ALPHA]]);
						gccp[BETA] = Double.parseDouble(cont[ind[BETA]]);
						gccp[OTHER] = Double.parseDouble(cont[ind[OTHER]]);
					}
				}
			}
			in.close();
			multiresult[0] = seq;
			multiresult[1] = go;
			multiresult[2] = ge;
			multiresult[3] = ssp;
			multiresult[4] = lccp;
			multiresult[5] = gccp;
		} catch (IOException e) {
			System.out.println("No Input (secondary structure preference)!");
		}
		return multiresult;
	}
	
	//potential reader
	public static double[][][] readPotentials(String helixFilepath, String betaFilepath, String coilFilepath){
		double[][][] potentials;
		try {
			//readin ccpa
			BufferedReader in = new BufferedReader(new FileReader(helixFilepath));
			String line;
			line = in.readLine();
			String temp[] = line.split("\\s+");
			potentials = new double[3][temp.length-1][26];
			while((line = in.readLine()) != null){
				temp = line.split("\\s+");
				char c = Bio.codeTranslate(temp[0]);
				for (int i = 0; i < potentials[0].length; i++) {
					potentials[0][i][c-65] = Double.parseDouble(temp[i+1]);
				}
			}
			in.close();
			
			//readin ccpb
			in = new BufferedReader(new FileReader(betaFilepath));
			line = in.readLine();
			while((line = in.readLine()) != null){
				temp = line.split("\\s+");
				char c = Bio.codeTranslate(temp[0]);
				for (int i = 0; i < potentials[0].length; i++) {
					potentials[1][i][c-65] = Double.parseDouble(temp[i+1]);
				}
			}
			in.close();
			
			//readin ccpo
			in = new BufferedReader(new FileReader(coilFilepath));
			line = in.readLine();
			while((line = in.readLine()) != null){
				temp = line.split("\\s+");
				char c = Bio.codeTranslate(temp[0]);
				for (int i = 0; i < potentials[0].length; i++) {
					potentials[2][i][c-65] = Double.parseDouble(temp[i+1]);
				}
			}
			in.close();
			return potentials;
		} catch (IOException e) {
			System.out.println("Cannot read potential file");
		}
		return null;
	}
	

	private static int[] calcIndices(String[] split) {
		int[] result = new int[3];
		for (int i = 1; i < split.length; i++) {
			if (split[i].equals("l"))
				result[OTHER] = i - 1;
			if (split[i].equals("h"))
				result[ALPHA] = i - 1;
			if (split[i].equals("s"))
				result[BETA] = i - 1;
		}
		return result;
	}

}
