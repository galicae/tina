package bioinfo.proteins.fr4gment;

import java.io.BufferedReader;
import java.io.FileReader;
import java.util.HashMap;
import java.util.LinkedList;

import bioinfo.proteins.PDBEntry;
import bioinfo.proteins.PDBFileReader;

/**
 * this class reads a cluster file containing a multiple superposition. Credit
 * for the construction of the set goes to Marie and the girl group. Credit for
 * not constructing it in a friendly way goes to them as well. (edit seitza):
 * TRUE STORY! =)
 * 
 * @author galicae
 * 
 */
public class MultipleSuperposition {
	private LinkedList<PDBEntry> models;
	private String directory;
	private int size;

	public MultipleSuperposition(String directory) {
		models = new LinkedList<PDBEntry>();
		this.directory = directory;
		readMultipleSuperpos();
		this.size = models.size();
	}

	public int length() {
		return this.size;
	}

	public void setSize() {
		this.size = models.size();
	}

	/**
	 * filters for all structures from the same family
	 * 
	 * @param guideFile
	 * @param queryId
	 */
	public void filterFamily(String guideFile, String queryId) {
		HashMap<String, String> idMap = new HashMap<String, String>();
		try {
			BufferedReader r = new BufferedReader(new FileReader(guideFile));
			String line = "";

			while ((line = r.readLine()) != null) {
				for (PDBEntry e : models) {
					if (line.startsWith(e.getId())) {
						idMap.put(e.getId(), line.split("\\s+")[1]);
						continue;
					}
					if (line.startsWith(queryId)) {
						idMap.put(queryId, line.split("\\s+")[1]);
					}
				}
			}
			r.close();

			String queryFamily = idMap.get(queryId);
			for (int i = 0; i < models.size();) {
				PDBEntry e = models.get(i);
				if (queryFamily.equals(idMap.get(e.getId()))) {
					models.remove(e);
					continue;
				}
				i++;
			}
		} catch (Exception e) {
			e.printStackTrace();
		}

	}

	/**
	 * 
	 * @param guideFile
	 * @param queryId
	 */
	public void filterSuperfamily(String guideFile, String queryId) {
		HashMap<String, String> idMap = new HashMap<String, String>();
		try {
			BufferedReader r = new BufferedReader(new FileReader(guideFile));
			String line = "";

			while ((line = r.readLine()) != null) {
				for (PDBEntry e : models) {
					if (line.startsWith(e.getId())) {
						String cath = line.split("\\s+")[1];
						cath = getSuperfamily(cath);
						idMap.put(e.getId(), cath);
						break;
					}
					if (line.startsWith(queryId)) {
						String cath = line.split("\\s+")[1];
						cath = getSuperfamily(cath);
						idMap.put(queryId, cath);
						break;
					}
				}
			}
			r.close();

			String queryFamily = idMap.get(queryId);
			for (int i = 0; i < models.size();) {
				PDBEntry e = models.get(i);
				if (queryFamily.equals(idMap.get(e.getId()))) {
					models.remove(e);
					continue;
				}
				i++;
			}
		} catch (Exception e) {
			e.printStackTrace();
		}
	}

	/**
	 * parse a cath id to read the superfamily
	 * 
	 * @param cath
	 *            the cath id
	 * @return the superfamily (the whole id apart from the part after the last
	 *         dot)
	 */
	private String getSuperfamily(String cath) {
		String[] cathAn = cath.split("\\.");
		cath = cathAn[0] + "." + cathAn[1] + "." + cathAn[2];
		return cath;
	}

	/**
	 * this function reads the cluster file (a PDB file with multiple models)
	 * specified in {@value directory} and parses its contents.
	 */
	private void readMultipleSuperpos() {
		try {
			// first read the file and save each model as a string
			BufferedReader reader = new BufferedReader(
					new FileReader(directory));
			PDBFileReader pdbReader = new PDBFileReader();
			String line = "";
			String nextID = null;
			while ((line = reader.readLine()) != null) {
				if (line.startsWith("ENDMDL")) {
					continue;
				}
				if (line.startsWith("REMARK")) {
					nextID = line.split("\\s+")[1];
				}
				if (line.startsWith("MODEL")) {
					// and while reading also parse PDBs
					PDBEntry e = pdbReader.readPDBFromModel(nextID, reader);
					models.add(e);
				}
			}
			reader.close();
		} catch (Exception e) {
			e.printStackTrace();
			// System.out.println("couldn't find " + directory);
		}
	}

	public LinkedList<PDBEntry> getStructures() {
		return this.models;
	}

	public void sort(String id) {
		int index = -1;
		for (PDBEntry p : models) {
			index = models.indexOf(p);
			if (id.contains(p.getId())) {
				models.add(0, p);
				break;
			}
		}
		models.remove(index + 1);
	}
}

// MODEL 1m6bA03
// MODEL 1m6bA01
// MODEL 1kooB00
// MODEL 1ds9A00
// MODEL 1jl5A00
// MODEL 1ogqA00
// MODEL 1k5dC00
// MODEL 1xecB00
// MODEL 1p8tA00
// MODEL 1otoA00
// MODEL 1dceC03
// MODEL 1w8aA00
// MODEL 1ookG00
// MODEL 1fs2C00
// MODEL 1pgvA00
// MODEL 1io0A00
// MODEL 1a9nC00
// MODEL 1a4yD00

