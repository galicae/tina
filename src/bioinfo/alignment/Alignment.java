package bioinfo.alignment;

public interface Alignment {
	
	/**
	 * Calculates an 2D array similiar to the map[][], but that only contains
	 * the aligned residues
	 * @return the 2D array
	 */
	int[][] getAlignedResidues();
	/**
	 * calculates the map between seq1 and seq2.<br />
	 * map[i][j] is -1 if j-th character of row(i) is aligned to a gap and the
	 * number of the character of row(1-i) it is aligned to else.
	 * 
	 * @return the map between seq1 and seq2
	 */
	int[][] calcMap();
	/**
	 * The length of the Alignment
	 */
	int length();
	
	double getScore();
	String toString();
	String toStringVerbose();
	Alignable getComponent(int number);
	Alignment duplicate();
	
}
