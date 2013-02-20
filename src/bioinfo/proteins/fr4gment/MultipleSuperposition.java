package bioinfo.proteins.fr4gment;

import java.io.BufferedReader;
import java.io.FileReader;
import java.util.LinkedList;

import bioinfo.proteins.PDBEntry;
import bioinfo.proteins.PDBFileReader;

/**
 * this class reads a cluster file containing a multiple superposition. Credit
 * for the construction of the set goes to Marie and the girl group. Credit for
 * not constructing it in a friendly way goes to them as well.
 * 
 * @author galicae
 * 
 */
public class MultipleSuperposition {
	private LinkedList<PDBEntry> models;
	private String directory;

	public MultipleSuperposition(String directory) {
		models = new LinkedList<PDBEntry>();
		this.directory = directory;

		readMultipleSuperpos();
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
			while ((line = reader.readLine()) != null) {
				if (line.startsWith("ENDMDL")) {
					continue;
				}
				if (line.startsWith("MODEL")) {
					// and while reading also parse PDBs
					PDBEntry e = pdbReader.readPDBFromModel(
							line.split("\\s+")[1], reader);
					models.add(e);
				}
			}
			reader.close();
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
	
	public LinkedList<PDBEntry> getStructures() {
		return this.models;
	}
	
	public void sort(String id) {
		int index = -1;
		for(PDBEntry p: models) {
			index = models.indexOf(p);
			if(id.contains(p.getID())) {
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