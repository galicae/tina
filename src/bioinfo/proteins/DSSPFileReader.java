package bioinfo.proteins;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.List;

/**
 * DSSPFileReader read DSSPFile and returns it as internal DSSPEntry
 * 
 * @author andreseitz
 * 
 */
public class DSSPFileReader {

	/**
	 * here you can find the DSSP Files
	 */
	public static final String DSSP_FOLDER = "/home/h/huberste/gobi/data/dssp/";

	private String folder = null;
	private List<String> files = null;
	private int sequentialReadIndex = 0;

	/**
	 * standard constructor creating "empty" DSSPFileReader
	 */
	public DSSPFileReader() {

	}

	/**
	 * checked folder will end with "/"
	 * 
	 * @param folder
	 *            which ends with "/"
	 */
	public DSSPFileReader(String folder) {
		this.folder = checkfolder(folder);
	}

	/**
	 * 
	 * @param filename
	 *            representing DSSP file to read in format: relpath/aaaaA00.DSSP
	 * @return
	 */
	public DSSPEntry readDSSPFromFile(String filename) {
		BufferedReader br = null;
		try {
			br = new BufferedReader(new InputStreamReader(new FileInputStream(
					filename)));
		} catch (Exception e) {
			e.printStackTrace();
		}
		String id = filename.split("/")[filename.split("/").length - 1]
				.split("\\.")[0];
		return parseEntry(br, id);
	}

	/**
	 * @return List of DSSPEntries in folder, null if no DSSPs are in folder or
	 *         folder is not set
	 */
	public List<DSSPEntry> readDSSPFolder() {
		List<DSSPEntry> result = new ArrayList<DSSPEntry>();
		BufferedReader br = null;
		String id;
		try {
			if (folder == null) {
				return null;
			}
			File[] files = new File(folder).listFiles();
			if (files.length == 0) {
				return null;
			}
			for (File file : files) {
				br = new BufferedReader(new InputStreamReader(
						new FileInputStream(file)));
				id = file.getName().split("\\.")[0];
				result.add(parseEntry(br, id));
			}
			return result;
		} catch (Exception e) {
			e.printStackTrace();
		}
		return null;
	}

	/**
	 * 
	 * @param folder
	 *            to read DSSPs from, folder will be set as
	 *            folder-class-variable
	 * @return List of DSSPEntries in folder, null if no DSSPs are in folder or
	 *         folder is not set
	 */
	public List<DSSPEntry> readDSSPFolder(String folder) {
		this.folder = checkfolder(folder);
		List<DSSPEntry> result = new ArrayList<DSSPEntry>();
		BufferedReader br = null;
		String id;
		try {
			if (folder == null) {
				return null;
			}
			File[] files = new File(folder).listFiles();
			if (files.length == 0) {
				return null;
			}
			for (File file : files) {
				br = new BufferedReader(new InputStreamReader(
						new FileInputStream(folder + file)));
				id = file.getName().split("\\.")[0];
				result.add(parseEntry(br, id));
			}
			return result;
		} catch (Exception e) {
			e.printStackTrace();
		}
		return null;
	}

	/**
	 * 
	 * @param id
	 *            of DSSP to read from
	 * @return DSSPentry representing the DSSP of given id, null if folder is
	 *         not set or certain DSSP does not exist
	 */
	public DSSPEntry readFromFolderById(String id) {
		if (folder == null) {
			return null;
		}
		// huberste: DSSPFiles normally are named without the ChainID!
		// String filename = folder+id.substring(0, 4).toUpperCase()+".dssp";
		// String filename = folder+id.substring(0, 4).toUpperCase()+".dssp";
		// seitza: but our f*cking files are ^^ searched the mistake in my code
		// a whole night
		// galicae: I had the same problem with my clusters - I wrote the
		// parseRealEntry for
		// this exact reason!
		String filename = folder + id + ".dssp";
		BufferedReader br = null;
		try {
			br = new BufferedReader(new InputStreamReader(new FileInputStream(
					folder + id + ".dssp")));
		} catch (Exception e) {
			return null;
		}
		return parseEntry(br, id);
	}

	/**
	 * initialise DSSPEntry list to read it out sequentially
	 */
	public void initSequentialFolderRead() {
		if (folder == null) {
			return;
		}
		File[] files = new File(folder).listFiles();
		List<String> names = new ArrayList<String>();
		for (File file : files) {
			names.add(file.getName());
		}
		this.files = names;
		sequentialReadIndex = 0;
	}

	/**
	 * 
	 * @return next DSSPEntry from initialisedList
	 */
	public DSSPEntry nextDSSP() {
		if (files.size() > sequentialReadIndex) {
			String file = files.get(sequentialReadIndex);
			BufferedReader br;
			try {
				br = new BufferedReader(new InputStreamReader(
						new FileInputStream(folder + file)));
				String id = file.split("\\.")[0];
				sequentialReadIndex++;
				return parseEntry(br, id);
			} catch (Exception e) {
				e.printStackTrace();
			}
			return null;
		} else {
			return null;
		}

	}

	/**
	 * 
	 * @return folder to read from
	 */
	public String getFolder() {
		return folder;
	}

	/**
	 * 
	 * @param folder
	 *            to read from
	 */
	public void setFolder(String folder) {
		this.folder = checkfolder(folder);
	}

	/**
	 * 
	 * @return list of files being set for read out
	 */
	public List<String> getFiles() {
		return files;
	}

	/**
	 * 
	 * @param files
	 *            which shall be read from folder if folder != null or from
	 *            anywhere in the file system if folder == null
	 * @return true if files were set correctly or false if setFiles failed
	 */
	public boolean setFiles(List<String> files) {
		if (folder == null) {
			return false;
		}
		File[] existing = new File(folder).listFiles();
		boolean found = false;
		for (String file : files) {
			found = false;
			for (int i = 0; i != existing.length; i++) {
				if (file.equals(existing[i].getName())) {
					found = true;
					break;
				}
			}
			if (!found) {
				return false;
			}
		}
		this.files = files;
		return true;
	}

	/**
	 * 
	 * @param folder
	 *            which should be checked for validity
	 * @return
	 */
	private String checkfolder(String folder) {
		if (folder.endsWith("/")) {
			return folder;
		} else {
			return folder + "/";
		}
	}

	/**
	 * 
	 * @param br
	 *            BufferedReader to read from
	 * @param DSSPId
	 *            DSSPId of format aaaaA00
	 * @return DSSPEntry representing the read file
	 */
	private static DSSPEntry parseEntry(BufferedReader br, String DSSPId) {
		String line;
		boolean sectionFlag = false;
		List<AminoAcidName> amino = new ArrayList<AminoAcidName>();
		List<Integer> resIndex = new ArrayList<Integer>();
		List<SecStructEight> secstr = new ArrayList<SecStructEight>();
		List<Integer> accessability = new ArrayList<Integer>();
		List<Double> phi = new ArrayList<Double>();
		List<Double> psi = new ArrayList<Double>();
		List<double[]> caTrace = new ArrayList<double[]>();
		double[] coord;

		try {

			while ((line = br.readLine()) != null) {
				if (sectionFlag) {
					if (line.charAt(13) == '!') {
						continue;
					}
					if (line.charAt(13) > 96) {
						amino.add(AminoAcidName.getAAFromOLC((char) (line
								.charAt(13) - 32)));
					} else {
						amino.add(AminoAcidName.getAAFromOLC(line.charAt(13)));
					}
					resIndex.add(Integer.parseInt(line.substring(5, 10).trim()));
					secstr.add(SecStructEight.getSSFromChar(line.charAt(16)));
					accessability.add(Integer.parseInt(line.substring(34, 38)
							.trim()));
					phi.add(Double.parseDouble(line.substring(103, 109).trim()));
					psi.add(Double.parseDouble(line.substring(109, 115).trim()));
					coord = new double[3];
					coord[0] = Double.parseDouble(line.substring(115, 122)
							.trim());
					coord[1] = Double.parseDouble(line.substring(122, 129)
							.trim());
					coord[2] = Double.parseDouble(line.substring(129, 136)
							.trim());
					caTrace.add(coord);
				}
				if (line.trim().startsWith("#")) {
					sectionFlag = true;
				}
			}
			br.close();
		} catch (Exception e) {
			e.printStackTrace();
		}
		return new DSSPEntry(DSSPId, amino, resIndex, secstr, accessability,
				phi, psi, caTrace);
	}

	/**
	 * Reads an DSSPEntry from a file without the DSSPFileReader around
	 * 
	 * @param filename
	 *            (format: xxxxA00.dssp)
	 * @return the DSSPEntry
	 */
	public static DSSPEntry readDSSPFile(String filename) {

		String[] tmp = filename.split("/");
		String temp = tmp[tmp.length-1];
		tmp = temp.split("\\.");
		String id = tmp[0];
		BufferedReader br = null;
		try {
			br = new BufferedReader(new InputStreamReader(new FileInputStream(
					filename)));
		} catch (IOException e) {
			System.out.println("Error reading DSSPFile: "+e.getLocalizedMessage());
			e.printStackTrace();
			return null;
		}
		return parseEntry(br, id);
	}

}
