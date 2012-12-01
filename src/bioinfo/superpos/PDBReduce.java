package bioinfo.superpos;

import bioinfo.alignment.Alignment;
import static bioinfo.proteins.AtomType.*;
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
			if (alignment.calcMap()[0][i] != -1) {
				try {
					alignmentCoordinates[0][j] = pdb1.getAminoAcid(i)
							.getAtomByType(CA).getPosition();
					alignmentCoordinates[1][j] = pdb2
							.getAminoAcid(alignment.calcMap()[0][i])
							.getAtomByType(CA).getPosition();
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
}
