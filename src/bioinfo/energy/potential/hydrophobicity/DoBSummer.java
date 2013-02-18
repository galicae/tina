/******************************************************************************
 * bioinfo.energy.potential.hydrophobicity.DoBSummer.java                     *
 *                                                                            *
 * This class's main method adds the frequency of specified dob (degree of    *
 * burial) intervals for all  AminoAcidTypes over a list of DoBFreqCounter    *
 * output files.                                                              *
 *                                                                            *
 * This file is best read at line width 80 and tab width 4.                   *
 *                                                                   huberste *
 ******************************************************************************/
package bioinfo.energy.potential.hydrophobicity;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;

/**
 * 
 * @author huberste
 * @lastchange 2013-02-18
 */
public class DoBSummer {
	
	public final static String usage = 
		"usage:\n" +
		"\tjava DoBFreqCounter <infolder>";
	
	private final static int BREAKS = 1024;
	
	/**
	 * @param args infolder
	 */
	public static void main(String[] args) {
		if (args.length < 1) {
			System.out.println(usage);
			System.exit(1);
		}
		String infolder = args[0];
		
		BufferedReader br = null;
		File dir = null;
		long[][] freq = new long[26][BREAKS];
		try {
		
			dir = new File(infolder);
			String[] names = dir.list();
	
			// At this point Java listed all files in c:/tmp and loaded their names into an array of Strings
			for (String filename : names) {
				System.out.println(filename);
				String line = null;
				try {
					br = new BufferedReader(new FileReader(infolder+filename));
					int aa = 0;
					while ((line = br.readLine()) != null) {
						if (line.startsWith("#")) { // comment
							continue;
						}
						if (!line.isEmpty()) {	// e.g. last line
							aa = (int) line.charAt(0) -65;
							String[] temp = line.split("\t");
							for(int i = 0; i < BREAKS; i++) {
								freq[aa][i] += Long.parseLong(temp[i+1]);
							}
						}
					}
				} catch (IOException e) {
					System.err.println("Error 93: problems reading the File:");
					e.printStackTrace();
				} finally {
					try {
						if (br != null) {
							br.close();
							br = null; // give br to GarbageCollector
						}
					} catch (IOException e){
						System.err.println("Error 102: problems closing the File:");
						e.printStackTrace();
					}
				}
			}
		} catch (NullPointerException e) {
			
		} finally {
			try {
				if (dir != null)
					dir = null;
			} catch (Exception e) {
				System.err.println("Error" + e.getMessage() +":");
				e.printStackTrace();
			}
		}
		
		// print out frequencies
		char x = 0;
		for (int i = 0; i < 26; i++) {
			x = (char)(65+i);
			System.out.print(x);
			for (int j = 0; j < BREAKS; j++) {
				System.out.print("\t" + freq[i][j]);
			}
			System.out.print("\n");
		}

	}

}

/******************************************************************************
 * "A question that sometimes drives me hazy:                                 *
 *  Am I or are the others crazy?"                                            *
 *     - Albert Einstein (1879 - 1955)                                        *
 ******************************************************************************/
