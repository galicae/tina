/******************************************************************************
 * huberdp.parameters.MaxFinder.java                                          *
 *                                                                            *
 * Contains a main method to find the maximum scores of the output of         *
 * RDPParameterTester.                                                        *
 *                                                                            *
 * This file is best read at line width 80 and tab width 4.                   *
 *                                                                   huberste *
 ******************************************************************************/
package huberdp.parameters;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;

/**
 * @author huberste
 * @lastchange 2013-02-22
 */
public class MaxFinder {

	private static final String outputfile = "/home/h/huberste/gobi/out/out";
	private static final int MAX_NUM = 10;
	
	/**
	 * @param args
	 */
	public static void main(String[] args) {
		
		String[] maxLines = new String[MAX_NUM];
		double[] maxVals = new double[MAX_NUM];
		for (int i = 0; i < MAX_NUM; i++) maxVals[i] = 0;
		int min = 0;
		
		BufferedReader br = null;
		try {
			br = new BufferedReader(new FileReader(outputfile));
			String line = null;
			while ((line = br.readLine()) != null) {
				String[] temp = line.split(" ");
				if (Double.parseDouble(temp[8]) > maxVals[min]) {
					maxLines[min] = line;
					maxVals[min] = Double.parseDouble(temp[7]);
					for (int i = 0; i < MAX_NUM; i++) {
						if (maxVals[i] < maxVals[min]) {
							min = i;
						}
					}
				}
				
			}
			for (int i = 0; i < MAX_NUM; i++) {
				System.out.println(maxLines[i]);
			}
		} catch (IOException e) {
			System.err.println("Error reading file: "+e.getLocalizedMessage());
			e.printStackTrace();
		} finally {
			try {
				if (br != null) {
					br.close();
					br = null; // GC
				}
			} catch (IOException e) {
				System.err.println("Error closing file: "+e.getLocalizedMessage());
				e.printStackTrace();
			}
		}

	}

}
