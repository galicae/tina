/**
 * 
 */
package bioinfo.pdb;

import java.io.BufferedInputStream;
import java.io.BufferedOutputStream;
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.net.URL;

/**
 * @author huberste
 *
 */
public class PDBFile {

	public final static String PDB_FILE_URL = "http://www.rcsb.org/pdb/files/";
	public final static String PDB_FILE_ENDING = ".pdb";
	
//	private String filePath;
//	
//	public PDBFile(String pdbFilePath, String pdbID) {
		
//	}
	
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
	
	public static String downloadPDB (String pdbID, String toPath) {
		
		BufferedInputStream in = null;
		BufferedOutputStream out = null;
		try {
			URL url = new URL(PDB_FILE_URL+pdbID.toUpperCase()+PDB_FILE_ENDING);
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
