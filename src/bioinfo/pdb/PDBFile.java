/******************************************************************************
 * bioinfo.proteins.PDBFile.java                                              *
 * Provides some useful methods for PDB Files.                                *
 *                                                                            *
 * This file is best read at line width 80 and tab width 4.                   *
 *                                                                   huberste *
 ******************************************************************************/
package bioinfo.pdb;

import java.io.BufferedInputStream;
import java.io.BufferedOutputStream;
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.net.URL;

/**
 * @author huberste
 * @lastchange 2013-02-17
 */
public class PDBFile {

	public final static String PDB_FILE_URL = "http://www.rcsb.org/pdb/files/";
	public final static String PDB_FILE_ENDING = ".pdb";
	
	/**
	 * if the given PDB file exists in pdbFilePath this method just returns the
	 * file path, else it downloads the file and returns the downloaded file's
	 * path 
	 * @param pdbFilePath
	 * @param pdbID
	 * @return the PDB file's path that contains the requested PDB
	 */
	public static String getFile(String pdbFilePath, String pdbID) {
		String result = null;
		String filePath = pdbFilePath+pdbID.toUpperCase()+PDB_FILE_ENDING;
		if (new File(filePath).exists()) {
			result = filePath;
		} else {
			// Download PDBFile
			result = downloadPDB(pdbID, pdbFilePath);
		}
		
		return result;
	}
	
	/**
	 * Downloads a given PDB File
	 * @param pdbID PDB ID of the requested PDB
	 * @param toPath Path where the file shall be downloaded to
	 * @return the Path to the downloaded file
	 */
	public static String downloadPDB (String pdbID, String toPath) {
		
		BufferedInputStream in = null;
		BufferedOutputStream out = null;
		try {
			//TODO: remember to stash
			URL url = new URL(PDB_FILE_URL+pdbID.substring(0,4).toUpperCase()+PDB_FILE_ENDING);
			in = new BufferedInputStream(url.openStream());
			out = new BufferedOutputStream(new FileOutputStream(toPath+pdbID.toUpperCase()+PDB_FILE_ENDING));
			byte[] buff = new byte[1024];
			int count = 0;
			while((count = in.read(buff, 0, 1024)) > 0) {
				out.write(buff, 0, count);
			}
			
		} catch (Exception e) {
			// MalformedURLException
			// IOException
			e.printStackTrace();
		} finally {
			try {
				in.close();
				out.close();
			} catch (IOException e) {
				System.err.println("Error closing the Streams.");
				e.printStackTrace();
			}
		}
		
		return toPath+pdbID.toUpperCase()+PDB_FILE_ENDING;
		
	}
	
}

/******************************************************************************
 * "A question that sometimes drives me hazy:                                 *
 *  Am I or are the others crazy?"                                            *
 *     - Albert Einstein (1879 - 1955)                                        *
 ******************************************************************************/
