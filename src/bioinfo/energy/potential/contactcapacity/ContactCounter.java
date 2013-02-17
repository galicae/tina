/******************************************************************************
 * bioinfo.energy.potential.contactcapacity.ContactCounter.java               *
 *                                                                            *
 * This class's main method counts the (local and global) contacts for all    *
 * AminoAcidTypes over a list of pdb files.                                   *
 *                                                                            *
 * This file is best read at line width 80 and tab width 4.                   *
 *                                                                   huberste *
 ******************************************************************************/
package bioinfo.energy.potential.contactcapacity;

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
 * @lastchange 2013-02-17
 */
public class ContactCounter {
	
	public final static String usage = 
		"usage:\n" +
		"\tjava ContactCounter <pdblist> <pdbpath>\n\n" +
		"where <pdblist> is a list of PDB IDs, <pdbpath> is a\n"+
		"(writable) path to the directory containing the PDB files";
	
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
				if (!line.isEmpty())	// e.g. last line
					pdbIDs.add(line);
			}
		} catch (IOException e) {
			System.err.println("Error 65: problems reading the File:");
			e.printStackTrace();
		} finally {
			try {
				if (br != null) {
					br.close();
					br = null; // give br to GarbageCollector
				}
			} catch (IOException e){
				System.err.println("Error 74: problems closing the File:");
				e.printStackTrace();
			}
		}
		
		PDBFileReader pdbreader = new PDBFileReader(pdbpath);
		
		// TODO count the occurences of the AminoAcids with some specific
		// contact counts
//		long[][] freq = new long[26][breaks];
		
		int debug = 0;
		
		for (String id : pdbIDs) { // for each file
			debug++;
			try {
				// begin debugging
				System.out.println("working on id "+id + " ("+debug+" of " + pdbIDs.size()+")");
				// end debugging
				
				PDBEntry structure = pdbreader.readPDBFromFile(
						PDBFile.getFile(pdbpath, id.substring(0,4)),id.charAt(4));
				
				if (structure == null) { // dirty debugging
					System.err.println("Error 98: Error getting ID "+id);
					continue;
				}
				
				// for each amino acid in that pdb
				for (int pos = 0; pos < structure.length(); pos++) {
					int astype = (structure.getAminoAcid(pos).getName().getOneLetterCode().charAt(0))-65;
					// calculate number of contacts
					int local = 0;
					int longrange = 0;
				} // for each amino acid in that pdb
			} catch (Exception e) { // dirty debugging
				System.err.println("problem at id "+id + " ("+debug+" of " + pdbIDs.size()+")");
				e.printStackTrace();
			}
		} // for each pdb
		
		// print out frequencies
		char x = 0;
		for (int i = 0; i < 26; i++) {
			x = (char)(65+i);
			System.out.print(x);
//			for (int j = 0; j < breaks; j++) {
//				System.out.print("\t" + freq[i][j]);
//			}
			System.out.print("\n");
		}

	}

}

/******************************************************************************
 * "A question that sometimes drives me hazy:                                 *
 *  Am I or are the others crazy?"                                            *
 *     - Albert Einstein (1879 - 1955)                                        *
 ******************************************************************************/
