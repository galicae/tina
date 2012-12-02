package bioinfo.superpos;

import static bioinfo.proteins.AtomType.CA;

import java.util.LinkedList;

import bioinfo.alignment.Alignment;
import bioinfo.proteins.AminoAcid;
import bioinfo.proteins.Atom;
import bioinfo.proteins.PDBEntry;

/**
 * 
 * @author gobi4 PDBReduce contains only static reducing method which reduces
 *         two PDBs to double[][][] using an sequence-alignment
 * 
 */
public class PDBReduce {

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
		for (int i = 0; i < alignment.length(); i++) {
			if (calcMap[0][i] != -1) {
				alignmentLength++;
			}
		}
		// now go through the alignment again and copy the coordinates
		double[][][] alignmentCoordinates = new double[2][alignmentLength][3];
		int j = 0; // the index for the aligned positions
		for (int i = 0; i < alignment.length(); i++) {
			if (calcMap[0][i] != -1) {
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

	public static void reducePDB(Alignment alignment, PDBEntry pdb1, PDBEntry pdb2) {
		// first find out how many residues are aligned
		int alignmentLength = 0;
		int[][] calcMap = alignment.calcMap();
		for (int i = 0; i < alignment.length(); i++) {
			if (calcMap[0][i] != -1) {
				alignmentLength++;
			}
		}
		LinkedList<AminoAcid> pdb1List = new LinkedList<AminoAcid>();
		LinkedList<AminoAcid> pdb2List = new LinkedList<AminoAcid>();

		// now go through the alignment again and find all aligned amino acids
		// with a CA atom

		Atom[] temp = new Atom[1];
		for (int i = 0; i < alignment.length(); i++) {
			if (calcMap[0][i] != -1) {
				try {
					temp[0] = pdb1.getAminoAcid(i).getAtomByType(CA);
					AminoAcid aa1 = new AminoAcid(pdb1.getAminoAcid(i)
							.getName(), temp);
					pdb1List.addLast(aa1);
					temp[0] = pdb2.getAminoAcid(calcMap[0][i])
							.getAtomByType(CA);
					AminoAcid aa2 = new AminoAcid(pdb2.getAminoAcid(
							calcMap[0][i]).getName(), temp);
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
		pdb1 = new PDBEntry(pdb1.getID(), pdb1Array);
		pdb2 = new PDBEntry(pdb2.getID(), pdb2Array);

	}
}
