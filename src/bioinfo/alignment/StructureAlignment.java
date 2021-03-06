package bioinfo.alignment;

import java.text.DecimalFormat;
import java.util.List;

import bioinfo.proteins.AminoAcid;
import bioinfo.proteins.PDBEntry;

/**
 * 
 * @author gobi_4
 * @date 28.11.2012
 * 
 */
public class StructureAlignment implements Alignment {

	private static final DecimalFormat DFLONG = new DecimalFormat("0.0000");
	private static final DecimalFormat DFSHORT = new DecimalFormat("0.000");

	/**
	 * Sequence one; the top Sequence
	 */
	private final PDBEntry seq1;
	/**
	 * Sequence two; the bottom Sequence
	 */
	private final PDBEntry seq2;
	/**
	 * The two rows of the Alignment
	 */
	private final AminoAcid[][] rows;
	/**
	 * The score of the Alignment
	 */
	private final double score;
	/**
	 * The length of the Alignment
	 */
	private final int length;
	/**
	 * The indices of the aligned residues, used by Paul format: map is an array
	 * of arrays whereas the first dimension just differs between sequence 1 or
	 * 2 the scd dimension is an array of length n or m referring to the length
	 * n or m of the sequences 1 or 2. every position in the scd dimension of the
	 * arrays contains the index of the residue from the other sequence it is
	 * aligned to
	 * 
	 */
	private final int[][] map;

	/**
	 * Constructs an Alignment of seq1 and seq2 with the given map and score.
	 * 
	 * @param seq1
	 *            the 1st Sequence
	 * @param seq2
	 *            the 2nd Sequence
	 * @param rows
	 *            the two rows of the alignment
	 * @param score
	 *            the score of the Alignment
	 */
	public StructureAlignment(PDBEntry seq1, PDBEntry seq2, AminoAcid[][] rows,
			double score) {
		this.seq1 = seq1;
		this.seq2 = seq2;
		this.rows = rows;
		this.score = score;
		this.length = calcLength();
		this.map = calcMap();
	}

	/**
	 * Constructs an Alignment of seq1 and seq2 with the given alignment rows
	 * and score.
	 * 
	 * @param seq1
	 *            the 1st Sequence
	 * @param seq2
	 *            the 2nd Sequence
	 * @param row1
	 *            the first row of the Alignment
	 * @param row2
	 *            the second row of the Alignment
	 * @param score
	 *            the score of the Alignment
	 */

	public StructureAlignment(PDBEntry seq1, PDBEntry seq2, AminoAcid[] row1,
			AminoAcid[] row2, double score) {
		this.seq1 = seq1;
		this.seq2 = seq2;
		this.score = score;
		AminoAcid[][] temp = new AminoAcid[2][];
		temp[0] = row1;
		temp[1] = row2;
		this.rows = temp;
		this.length = row1.length;
		this.map = calcMap();
	}

	/**
	 * Constructs an Alignment of seq1 and seq2 with the given alignment rows
	 * and score.
	 * 
	 * @param seq1
	 *            the 1st Sequence
	 * @param seq2
	 *            the 2nd Sequence
	 * @param row1
	 *            the first row of the Alignment
	 * @param row2
	 *            the second row of the Alignment
	 * @param score
	 *            the score of the Alignment
	 */
	public StructureAlignment(PDBEntry seq1, PDBEntry seq2,
			List<AminoAcid> row1, List<AminoAcid> row2, double score) {
		this.seq1 = seq1;
		this.seq2 = seq2;
		this.score = score;
		AminoAcid[][] temp = new AminoAcid[2][];
		temp[0] = row1.toArray(new AminoAcid[row1.size()]);
		temp[1] = row2.toArray(new AminoAcid[row2.size()]);
		this.rows = temp;
		this.length = row1.size();
		this.map = calcMap();
	}

	/**
	 * Calculates the length of the given alignment.
	 * 
	 * @return the length of the alignment
	 */
	private int calcLength() {
		return rows[0].length;
	}

	/**
	 * calculates the map between seq1 and seq2.<br />
	 * map[i][j] is -1 if j-th character of row(i) is aligned to a gap and the
	 * number of the character of row(1-i) it is aligned to else.
	 * 
	 * @return the map between seq1 and seq2
	 */
	public int[][] calcMap() {
		int[][] result = new int[2][];
		result[0] = new int[seq1.length()];
		result[1] = new int[seq2.length()];
		int col1 = 0;
		int col2 = 0;
		for (int i = 0; i < rows[0].length; i++) {
			if (rows[0][i] == null) { // insertion
				result[1][col2] = -1;
				col2++;
			} else if (rows[1][i] == null) { // deletion
				result[0][col1] = -1;
				col1++;
			} else { // (mis)match
				result[0][col1] = col2;
				result[1][col2] = col1;
				col1++;
				col2++;
			}
		}
		return result;
	}

	/**
	 * 
	 * @return String of the format "idone idtwo score"
	 */
	public String toString() {
		String result = seq1.getId() + " " + seq2.getId() + " "
				+ DFLONG.format(score);
		return result;
	}

	/**
	 * 
	 * @return String of the format ">idone idtwo score\n idone: rowone\n idtwo:
	 *         rowtwo"
	 */
	public String toStringVerbose() {
		String amino1 = "";
		String amino2 = "";

		for (int i = 0; i != rows[0].length; i++) {
			if (rows[0][i] == null) {
				amino1 += "-";
			} else {
				amino1 += rows[0][i].toString();
			}
			if (rows[1][i] == null) {
				amino2 += "-";
			} else {
				amino2 += rows[1][i].toString();
			}
		}

		String result = ">" + seq1.getId() + " " + seq2.getId() + " "
				+ DFSHORT.format(score) + "\n" + seq1.getId() + ": " + amino1
				+ "\n" + seq2.getId() + ": " + amino2;
		return result;
	}

	/**
	 * 
	 * @return length of the Alignment
	 */
	public int length() {
		return length;
	}

	public double getScore() {
		return score;
	}

	public PDBEntry getComponent(int n) {
		if (n == 0)
			return seq1;
		else
			return seq2;
	}

	public AminoAcid[] getRow(int n) {
		return rows[n];
	}

	public StructureAlignment duplicate() {
		return new StructureAlignment(seq1, seq2, rows, score);
	}

	@Override
	//get an array with indices of only alignedResidues
	public int[][] getAlignedResidues() {
		int[][] result = new int[2][countAlignedResidues()];
		int temp = 0;
		for (int i = 0; i < map[0].length; i++) {
			if(map[0][i] != -1){
				result[0][temp] = i;
				result[1][temp] = map[0][i];
				temp++;
			}
		}
		return result;
	}

	// calc aligned indices
	public int countAlignedResidues() {
		int result = 0;
		for (int i = 0; i < map[0].length; i++) {
			if (map[0][i] != -1) {
				result++;
			}
		}
		return result;
	}
}
