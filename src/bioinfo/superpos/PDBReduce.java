/******************************************************************************
 * bioinfo.superpos.RDBReduce.java                                            *
 *                                                                            *
 * This file is best read at line width 80 and tab width 4.                   *
 *                                                                      gobi4 *
 ******************************************************************************/
package bioinfo.superpos;

import static bioinfo.proteins.AtomType.CA;

import java.util.LinkedList;

import bioinfo.alignment.Alignment;
import bioinfo.proteins.AminoAcid;
import bioinfo.proteins.Atom;
import bioinfo.proteins.AtomType;
import bioinfo.proteins.PDBEntry;

/**
 * PDBReduce contains only static reducing functions which reduce PDBs. See each
 * function for details.
 * 
 * @author gobi4
 * @lastchange 2013-02-14
 */
public class PDBReduce {

	// TODO JavaDoc
	/**
	 * 
	 * @param alignment
	 * @param pdb1
	 * @param pdb2
	 * @return double[][][] containing two double[][] Ca coordinate matrices for
	 *         kabsch method
	 */
	public static double[][][] reduce(Alignment alignment, PDBEntry pdb1,
			PDBEntry pdb2) {
		// first find out how many residues are aligned
		int alignmentLength = 0;
		int[][] calcMap = alignment.calcMap();
		for (int i = 0; i < calcMap[0].length; i++) {
			if (calcMap[0][i] > 0) {
				alignmentLength++;
			}
		}
		// now go through the alignment again and copy the coordinates
		double[][][] alignmentCoordinates = new double[2][alignmentLength][3];
		int j = 0; // the index for the aligned positions
		for (int i = 0; i < calcMap[0].length; i++) {
			if (calcMap[0][i] >= 0) {
				try {
					alignmentCoordinates[0][j] = pdb1.getAminoAcid(i)
							.getAtomByType(CA).getPosition();
					alignmentCoordinates[1][j] = pdb2
							.getAminoAcid(calcMap[0][i]).getAtomByType(CA)
							.getPosition();
					j++;
				}
				// catch amino acids with no CA atom
				catch (Exception e) {
					continue;
				}
			}
		}
		return alignmentCoordinates;
	}

	/**
	 * alternative version of reducePDB, after we stopped reading files and
	 * started using objects instead
	 * 
	 * @param alignment
	 *            an Alignment object
	 * @param pdb1
	 *            the first pdb (p)
	 * @param pdb2
	 *            the second pdb (q)
	 * @return two reduced pdb entries, where only the aligned Ca atons
	 */
	public static PDBEntry[] reducePDB(Alignment alignment, PDBEntry pdb1,
			PDBEntry pdb2) {
		// first find out how many residues are aligned
//		int alignmentLength = 0;
		int[][] calcMap = alignment.calcMap();
//		for (int i = 0; i < calcMap[0].length; i++) {
//			if (calcMap[0][i] >= 0) {
//				alignmentLength++;
//			}
//		}
		LinkedList<AminoAcid> pdb1List = new LinkedList<AminoAcid>();
		LinkedList<AminoAcid> pdb2List = new LinkedList<AminoAcid>();

		// now go through the alignment again and find all aligned amino acids
		// with a CA atom

		Atom[] temp = new Atom[1];
		for (int i = 0; i < calcMap[0].length; i++) {
			if (calcMap[0][i] >= 0) {
				try {
					temp[0] = pdb1.getAminoAcid(i).getAtomByType(CA);
					AminoAcid aa1 = new AminoAcid(pdb1.getAminoAcid(i)
							.getName(), pdb2.getAminoAcid(calcMap[0][i])
							.getResIndex(), temp.clone());
					pdb1List.addLast(aa1);
					temp[0] = pdb2.getAminoAcid(calcMap[0][i])
							.getAtomByType(CA);
					AminoAcid aa2 = new AminoAcid(pdb2.getAminoAcid(
							calcMap[0][i]).getName(), pdb2.getAminoAcid(
							calcMap[0][i]).getResIndex(), temp.clone());
					pdb2List.addLast(aa2);
				} catch (Exception e) {
					continue;
				}
			}
		}
		// now make arrays out of the lists
		AminoAcid[] pdb1Array = new AminoAcid[pdb1List.size()];
		AminoAcid[] pdb2Array = new AminoAcid[pdb1List.size()];
		for (int i = 0; i < pdb1List.size(); i++) {
			pdb1Array[i] = pdb1List.get(i);
			pdb2Array[i] = pdb2List.get(i);
		}

		// now change the PDBEntries, so that we can write them in files again
		pdb1 = new PDBEntry(pdb1.getId(), pdb1Array);
		pdb2 = new PDBEntry(pdb2.getId(), pdb2Array);
		PDBEntry[] result = { pdb1, pdb2 };
		return result;
	}
	
	/**
	 * Reduces one PDBEntry to a double[][], containing only the coordinates
	 * of all CA atoms
	 * @param pdb1
	 * @return a double[n][x|y|z] containing the (x|y|z) coordinates of the n-th
	 * CA atom 
	 */
	public static double[][] reduceSinglePDB(PDBEntry pdb1) {
		double[][] result = new double[pdb1.length()][3];
		Atom tempA = new Atom(AtomType.CA, new double[3]);
		for(int i = 0; i < pdb1.length(); i++) {
			tempA = pdb1.getAminoAcid(i).getAtomByType(AtomType.CA);
			if(tempA == null)
				return null;
			result[i] = tempA.getPosition();
		}
		return result;
	}
	
	/**
	 * Reduces two PDBEntries of the same protein (e.g. two different models for
	 * one protein) by removing all points in each model that does not exist in
	 * the other one.
	 * @author huberste
	 * @param arg1
	 * @param arg2
	 * @return the reduced point groups with
	 * result[0] first PDBEntry, result[1] second PDBEntry,
	 * result[n][x] is the x-th atom of the n-th PDBEntry,
	 * result[n][x][x|y|z] is the coordinate value
	 * @throws Exception 
	 */
	public static double[][][] reducePDBs(PDBEntry arg1, PDBEntry arg2)
			throws Exception {
		
		// The following code is not necessary AND it produces errors!
//		if (arg1.length() != arg2.length() ) {
//			throw new Exception("Both given PDBEntries don't seem to be "+
//					"representing the same protein.");
//		}
		
		// result[0] first PDBEntry, result[1] second PDBEntry
		// result[n][x] is the x-th atom of the n-th PDBEntry
		// result[n][x][x|y|z] is the coordinate value
		double[][][] result = new double[2][][];
		LinkedList<double[]> atoms1 = new LinkedList<double[]>();
		LinkedList<double[]> atoms2 = new LinkedList<double[]>();
		
		// go through every AminoAcid
		for (int i = 0; i < arg1.length(); i++) {
			AminoAcid as1 = arg1.getAminoAcid(i);
			AminoAcid as2 = arg2.getAminoAcid(i);
			
			// go through every Atom
			for (int atom1 = 0; atom1 < as1.getAtomNumber(); atom1++) {
				for (int atom2 = 0; atom2 < as1.getAtomNumber(); atom2++) {
					// if same Atom
					if (as1.getAtom(atom1).getType() ==
						as2.getAtom(atom2).getType()) {
						// save Atom to AtomList
						atoms1.add(as1.getAtom(atom1).getPosition());
						atoms2.add(as2.getAtom(atom2).getPosition());
		}	}	}	}
		result[0] = atoms1.toArray(new double[0][]);
		result[1] = atoms2.toArray(new double[0][]);
		
		return result;
		
	}
	
}
