/******************************************************************************
 * files.PairFile                                                             *
 *                                                                            *
 * This file is best read at line width 80 and tab width 4.                   *
 *                                                                   huberste *
 ******************************************************************************/
package files;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.LinkedList;

/**
 * The class PairFile is for reading pairs of sequences (aka Joblist)
 * @author huberste
 * @date November 11, 2012
 * @version 1.0
 */
public class PairFile {
	
	private String filename;
	
	private LinkedList<String[]> joblist;
	
	/**
	 * constructs a new pairfile object. Doesn't load the file, this only
	 * happens in getJobList().
	 * @param filename
	 */
	public PairFile(String filename) {
		this.filename = filename;
	}
	
	/**
	 * reads the joblist from the file.
	 * @return the jobList als LinkedList.
	 */
	public LinkedList<String[]> getJoblist() {
		if (joblist == null) {
			joblist = new LinkedList<String[]>();
			// following lines only usable in java7
/*
			Path path = Paths.get(filename);
			Charset charset = Charset.forName("UTF-8");
*/
			BufferedReader reader = null;
			try {
				reader = new BufferedReader(new FileReader(filename));
				String line = null;
				while ((line = reader.readLine()) != null) {
					// TODO optimization
					String[] temp = line.split("[ \t]");
					String[] temp2 = new String[2];
					temp2[0] = temp[0]; temp2[1]=temp[1];
					joblist.add(temp2);
				}
			} catch (IOException x) {
				System.err.format("IOException: %s%n", x);
			} finally {
				try {
					reader.close();
				} catch (IOException x) {
					System.err.format("IOException: %s%n", x);
				}
			}
			
		}
		return joblist;
	}
	
}

/******************************************************************************
 * "A question that sometimes drives me hazy:                                 *
 *  Am I or are the others crazy?"                                            *
 *     - Albert Einstein (1879 - 1955)                                        *
 ******************************************************************************/