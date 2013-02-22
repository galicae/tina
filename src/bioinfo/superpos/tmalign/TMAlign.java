package bioinfo.superpos.tmalign;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;

public class TMAlign {

	/**
	 * (full) path to TMAlign binary
	 */
	String binpath;

	/**
	 * Constructs a new TMAlign Object
	 * 
	 * @param binpath
	 */
	public TMAlign(String binpath) {
		this.binpath = binpath;
	}

	/**
	 * 
	 * @param templatepdbfile
	 * @param targetpdbfile
	 * @return
	 */
	public String align(String templatepdbfile, String targetpdbfile) {

		String command = binpath + " " + templatepdbfile + " " + targetpdbfile;
		try {
			return bioinfo.proteins.fragm3nt.run.RunHelper
					.execToString(command);
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}

		return null;
	}

	/**
	 * 
	 * @param templatepdbfile
	 * @param targetpdbfile
	 * @param matrixfile
	 * @return
	 */
	public String align(String templatepdbfile, String targetpdbfile,
			String matrixfile) {

		String command = binpath + " " + templatepdbfile + " " + targetpdbfile
				+ " -m " + matrixfile;
		try {
			return bioinfo.proteins.fragm3nt.run.RunHelper
					.execToString(command);
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}

		return null;
	}

	public static double[][] loadRotationMatrixFromFile(String filename) {
		double[][] result = new double[4][3];
		BufferedReader br = null;
		try {
			br = new BufferedReader(new FileReader(filename));
			String line = null;
			br.readLine();	// skip 1st line
			br.readLine();	// skip 2nd line
			line = br.readLine();	// 3rd line contains stuff
			result[0][0] = Double.parseDouble(line.substring(3,20).trim());
			result[1][0] = Double.parseDouble(line.substring(21, 35).trim());
			result[2][0] = Double.parseDouble(line.substring(36, 50).trim());
			result[3][0] = Double.parseDouble(line.substring(51, 65).trim());
			line = br.readLine();	// 4rd line contains stuff
			result[0][1] = Double.parseDouble(line.substring(3, 20).trim());
			result[1][1] = Double.parseDouble(line.substring(21, 35).trim());
			result[2][1] = Double.parseDouble(line.substring(36, 50).trim());
			result[3][1] = Double.parseDouble(line.substring(51, 65).trim());
			line = br.readLine();	// 5rd line contains stuff
			result[0][2] = Double.parseDouble(line.substring(3, 20).trim());
			result[1][2] = Double.parseDouble(line.substring(21, 35).trim());
			result[2][2] = Double.parseDouble(line.substring(36, 50).trim());
			result[3][2] = Double.parseDouble(line.substring(51, 65).trim());
			
		} catch (IOException e) {
			System.err
					.println("Error reading file: " + e.getLocalizedMessage());
			e.printStackTrace();
		} finally {
			try {
				if (br != null) {
					br.close();
					br = null; // GC
				}
			} catch (IOException e) {
				System.err.println("Error closing file: "
						+ e.getLocalizedMessage());
				e.printStackTrace();
			}
		}
		return result;
	}

}
