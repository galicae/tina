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
import java.util.LinkedList;

import static bioinfo.energy.potential.contactcapacity.ContactCapacityMatrix.*;
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
	
	public final static String filedesc =
			"# Contact Count Matrices.\n" +
			"# to the right are the number of long range contacts,\n" +
			"# downwards the number of local contacts";
	
	public final static String LONGRANGE_HEADER =
			"#\t\t\t# of long range contacts ->";
	
	public final static String LOCAL_HEADER =
			"#\t| # of local\n"+
			"#\tV   contacts";
	
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
				if (line.trim().startsWith("#")) { // comment
					continue;
				}
				if (!line.trim().isEmpty())	// e.g. last line
					pdbIDs.add(line.trim());
			}
		} catch (IOException e) {
			System.err.println("Error 60: problems reading the File:");
			e.printStackTrace();
		} finally {
			try {
				if (br != null) {
					br.close();
					br = null; // give br to GarbageCollector
				}
			} catch (IOException e){
				System.err.println("Error 69: problems closing the File:");
				e.printStackTrace();
			}
		}
		
		PDBFileReader pdbreader = new PDBFileReader(pdbpath);
		
		// TODO count the occurences of the AminoAcids with some specific
		// contact counts
		// TODO Where to set the (virtual) C-beta?
		long[][][] counts = new long
				[26]
				[MAX_LOCAL_CONTACTS + 1]
				[MAX_LONGRANGE_CONTACTS + 1];
		
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
					for (int pos2 = 0; pos2 < structure.length(); pos2++) {
						if (pos != pos2) { // natürlich nicht bei selber AS
							// TODO Check if in contact distance
							
							
							if (Math.abs(pos-pos2) <= MAX_LOCAL_DISTANCE) { // local
								local++;
							} else { // longRange
								longrange++;
							}
						}
					}
					
					// insert into count matrix
					for (int i1 = MAX_LOCAL_CONTACTS; i1 >= 0; i1--) {
						if (local >= i1) {
							for (int i2 = MAX_LONGRANGE_CONTACTS; i2 >= 0; i2--) {
								if (longrange >= i2) {
									counts[astype][i1][i2]++;
									
									continue;
								}
							}
							
							continue;
						}
					}
					
				} // for each amino acid in that pdb
			} catch (Exception e) { // dirty debugging
				System.err.println("problem at id "+id + " ("+debug+" of " + pdbIDs.size()+")");
				e.printStackTrace();
			}
		} // for each pdb
		
		// print out frequencies
		
		System.out.println(filedesc);
		char x = 0;
		for (int aa = 0; aa < 26; aa++) {
			x = (char)(65+aa);
			System.out.print("# \n# " + x+ "\n# \n\n");
			System.out.println(LONGRANGE_HEADER);
			for (int local = 0; local <= MAX_LONGRANGE_CONTACTS; local++) {
				System.out.print("\t"+local);
			}
			System.out.print("+\n");
			
			for (int local = 0; local <= MAX_LOCAL_CONTACTS; local++) {
				System.out.print(local);
				for (int longrange = 0; longrange <= MAX_LONGRANGE_CONTACTS; longrange++) {
					System.out.print("\t" + counts[aa][local]);
				}
			}
			System.out.print("\n\n");
		}

	}
	
/* outputfile looks like this:
#
# A
# 

#			# of long range contacts ->
	0	1	2	3	4	5	6	7	8	9+
#	| # of local
#	V contacts
0	x	x	x	x	x	x	x	x	x	x
1	x	x	x	x	x	x	x	x	x	x
2	x	x	x	x	x	x	x	x	x	x
3	x	x	x	x	x	x	x	x	x	x
4	x	x	x	x	x	x	x	x	x	x
5	x	x	x	x	x	x	x	x	x	x
6	x	x	x	x	x	x	x	x	x	x
7	x	x	x	x	x	x	x	x	x	x
8	x	x	x	x	x	x	x	x	x	x
9	x	x	x	x	x	x	x	x	x	x
10	x	x	x	x	x	x	x	x	x	x
11	x	x	x	x	x	x	x	x	x	x
12+	x	x	x	x	x	x	x	x	x	x

...
 */

}

/******************************************************************************
 * "A question that sometimes drives me hazy:                                 *
 *  Am I or are the others crazy?"                                            *
 *     - Albert Einstein (1879 - 1955)                                        *
 ******************************************************************************/
