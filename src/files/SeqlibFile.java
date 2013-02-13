/******************************************************************************
 * files.SeqlibFile                                                           *
 *                                                                            *
 * This file is best read at line width 80 and tab width 4.                   *
 *                                                                   huberste *
 ******************************************************************************/
package files;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.HashMap;
import bioinfo.Sequence;

/**
 * The class SeqlibFile is for loading a Sequence Library from file into a
 * HashMap. The IDs of the Sequences will be the keys, the sequences itself the
 * values.
 * @author huberste
 * @date November 25, 2012
 * @version 1.0
 */
public class SeqlibFile {
	
	private String filename;
	private HashMap<String, Sequence> library;
	
	/**
	 * constructs a new SeqLibFile Object. Doesn't load the data from file, this
	 * only happens in getLibrary().
	 * @param filename
	 */
	public SeqlibFile(String filename) {
		this.filename = filename;
	}
	
	/**
	 * loads the library into a HashMap. The IDs of the Sequences will be the
	 * keys, the sequences itself the values.
	 * @return the HashMap with the sequences IDs as keys and the sequences
	 * as values
	 */
	public HashMap<String, Sequence> getLibrary() {
		// TODO optimization: maybe only read needed sequences?
		if (library == null) {
			library = new HashMap<String, Sequence>();
			// Folliwing lines only usable in Java7
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
					String[] temp= line.split(":");
					library.put(temp[0],new Sequence(temp[0], temp[1].toCharArray()));
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
		return library;
	}
	
}
