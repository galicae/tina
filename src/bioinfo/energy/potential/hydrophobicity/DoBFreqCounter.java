/******************************************************************************
 * bioinfo.energy.potential.hydrophobicity.DoBFreqCounter.java                *
 *                                                                            *
 * This class's main method counts the frequency of specified dob (degree of  *
 * burial) intervals for all  AminoAcidTypes over a list of pdb files.        *
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
 * @lastchange 2013-02-16
 */
public class DoBFreqCounter {
	
	public final static String usage = 
		"usage:\n" +
		"\tjava DoBFreqCounter <pdblist> <pdbpath> <breaks>\n\n" +
		"where <pdblist> is a list of PDB IDs, <pdbpath> is a\n"+
		"(writable) path to the directory containing the PDB files and" +
		"<breaks> is the number of breaks you wish to make (e.g. 1024)";
	
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
	 * empirically calibratet value for voro++
	 */
	private final static double MIN_CONTACT = 2.0;
	
	/**
	 * @param args pdblistpath, pdbpath
	 */
	public static void main(String[] args) {
		if (args.length < 3) {
			System.out.println(usage);
			System.exit(1);
		}
		String pdbList = args[0];
		String pdbpath = args[1];
		int breaks = Integer.parseInt(args[2]);
		
//		// initialize important stuff	
//		Locale.setDefault(Locale.US);
//		DecimalFormat df = new DecimalFormat("0.000000");
		
		// load file
		LinkedList<String> pdbIDs = new LinkedList<String>();
		BufferedReader br = null;
		String line = null;
		
		try {
			br = new BufferedReader(new FileReader(pdbList));
			while ((line = br.readLine()) != null) {
				if (line.startsWith("#")) { // comment
					continue;
				}
				pdbIDs.add(line);
			}
		} catch (IOException e) {
			System.err.println("Error 77: problems reading the File:");
			e.printStackTrace();
		} finally {
			try {
				if (br != null) {
					br.close();
					br = null; // give br to GarbageCollector
				}
			} catch (IOException e){
				System.err.println("Error 83: problems closing the File:");
				e.printStackTrace();
			}
		}
		
		PDBFileReader pdbreader = new PDBFileReader(pdbpath);
		
		// count the frequencies of the AminoAcids with some specific d.o.b.
		long[][] freq = new long[26][breaks];
		
		int debug = 0;
		
		for (String id : pdbIDs) { // for each file
			debug++;
			// begin debugging
			System.out.println("working on id "+id + " ("+debug+" of " + pdbIDs.size()+")");
			// end debugging
			
			PDBEntry structure = pdbreader.readPDBFromFile(
					PDBFile.getFile(pdbpath, id.substring(0,4)),id.charAt(4));
			
			if (structure == null) { // dirty debugging
				System.err.println("error getting ID "+id);
				continue;
			}
			
			VoroPPWrap voro = new VoroPPWrap(VOROPATH);
			VoronoiData data = new VoronoiData(structure.getID());
			data.reducePDB(VoroPrepType.CA, structure);
			data.fillGridWithoutClashes(GRID_EXTEND, GRID_DENSITY, GRID_CLASH);
			voro.decomposite(data);
			data.detectOuterGrid(MIN_CONTACT);
			Set <Integer> solvents = data.getOuterGridIds();
		
			// for each amino acid
			for (int pos = 0; pos < structure.length(); pos++) {
				int astype = (structure.getAminoAcid(pos).getName().getOneLetterCode().charAt(0))-65;
				// calculate dob (degree of burial)
				double outer = 0.0;
				double inner = 0.0;
				double dob;
				HashMap<Integer,Double> faces = data.getFaces().get(pos);
				
				if (faces != null) { // dirty debugging
				
					for(int neighbor : faces.keySet()){
						if(solvents.contains(neighbor)){
							outer += faces.get(neighbor);
						}else{
							inner += faces.get(neighbor);
						}
					}
					
					dob = outer/(outer+inner);
					
					for(int i = 0; i < breaks; i++) {
						if (dob <= ((double)i+1.0)*(1.0/(double)breaks)) {
							freq[astype][i] ++;
							break;
						}
					}
					
					
				}
			} // for each amino acid in that pdb 		
		} // for each pdb
		
		// print out frequencies
		char x = 0;
		for (int i = 0; i < 26; i++) {
			x = (char)(65+i);
			System.out.print(x);
			for (int j = 0; j < breaks; j++) {
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
