/******************************************************************************
 * huberdp.Scoring.RDPScoring.java                                            *
 * This file contains the class RDPScoring which is RDP's scoring function.   *
 *                                                                            *
 * This file is best read at line width 80 and tab width 4.                   *
 *                                                                   huberste *
 ******************************************************************************/
package bioinfo.energy.potential.hydrophobicity;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.Set;

import bioinfo.energy.potential.voronoi.VoroPPWrap;
import bioinfo.energy.potential.voronoi.VoroPrepType;
import bioinfo.energy.potential.voronoi.VoronoiData;
import bioinfo.pdb.PDBFile;
import bioinfo.proteins.PDBEntry;
import bioinfo.proteins.PDBFileReader;

/**
 * 
 * @author huberste
 * @lastchange 2013-02-15
 */
public class HydrophobicityCalculator {
	
	public final static String usage = 
		"usage:\n" +
		"\tjava HydrobhobicityCalculator <pdblist> <pdbpath>\n\n" +
		"where <pdblist> is a list of PDB IDs and <pdbpath> is a\n"+
		"(writable) path to the directory containing the PDB files";
	
	/**
	 * static reference to voro++ path
	 */
	private final static String VOROPATH =
			"/home/h/huberste/gobi/tina/tools/voro++_ubuntuquantal";
	/**
	 * empirically calibratet value for voro++
	 */
	private final static double GRID_EXTEND = 8.9;
	/**
	 * empirically calibratet value for voro++
	 */
	private final static double GRID_DENSITY = 1.0;
	/**
	 * empirically calibratet value for voro++
	 */
	private final static double GRID_CLASH = 6.5;
	
	/**
	 * @param args pdblistpath, pdbpath
	 */
	public static void main(String[] args) {
		if (args.length < 2) {
			System.out.println(usage);
			System.exit(1);
		}
		String pdbList = args[0];
		String pdbpath = args[1];
		
		// TODO load file pdb_25
		LinkedList<String> pdbIDs = new LinkedList<String>();
		BufferedReader br = null;
		String line = null;
		
		try {
			br = new BufferedReader(new FileReader(pdbList));
			while ((line = br.readLine()) != null) {
				pdbIDs.add(line);
			}
		} catch (IOException e) {
			System.err.println("Error 77: problems reading the File:");
			e.printStackTrace();
		} finally {
			try {
				br.close();
			} catch (IOException e){
				System.err.println("Error 83: problems closing the File:");
				e.printStackTrace();
			}
		}
		
		// for each file
		PDBFileReader pdbreader = new PDBFileReader(pdbpath);
		int[] count = new int[26];
		double[] hydro = new double[26];
		for (String id : pdbIDs) {
			PDBEntry structure = pdbreader.readPDBFromFile(
					PDBFile.getFile(pdbpath, id.substring(0,4)),id.charAt(4));
			
			VoroPPWrap voro = new VoroPPWrap(VOROPATH);
			VoronoiData data = new VoronoiData(structure.getID());
			data.reducePDB(VoroPrepType.CA, structure);
			data.fillGridWithoutClashes(GRID_EXTEND, GRID_DENSITY, GRID_CLASH);
			voro.decomposite(data);
			Set <Integer> solvents = data.getOuterGridIds();
		
			// for each amino acid
			for (int pos = 0; pos < structure.length(); pos++) {
				
				int astype = (structure.getAminoAcid(pos).getName().getOneLetterCode().charAt(0))-65;
				// calculate dob (degree of burial)
				double outer = 0.0;
				double inner = 0.0;
				HashMap<Integer,Double> faces = data.getFaces().get(pos);
				
				for(int neighbor : faces.keySet()){
					if(solvents.contains(neighbor)){
						outer += faces.get(neighbor);
					}else{
						inner += faces.get(neighbor);
					}
				}
				hydro[astype] += outer/(outer+inner);
				count[astype]++;
			} // for each amino acid in that pdb 		
		} // for each pdb
		
		for(int i = 0; i < 26; i++) {
			if (count[i] > 0)
				hydro[i] = hydro[i]/(double)count[i];
			System.out.println(hydro[i]);
		}
		
		// TODO Verify results
		

	}

}

/******************************************************************************
 * "A question that sometimes drives me hazy:                                 *
 *  Am I or are the others crazy?"                                            *
 *     - Albert Einstein (1879 - 1955)                                        *
 ******************************************************************************/
