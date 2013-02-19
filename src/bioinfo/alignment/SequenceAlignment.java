/******************************************************************************
 * bioinfo.alignment.Alignment                                                *
 *                                                                            *
 * This file is best read at line width 80 and tab width 4.                   *
 ******************************************************************************/
package bioinfo.alignment;

import java.text.DecimalFormat;

import bioinfo.Sequence;

/**
 * @author gobi_4
 * @lastchange 2013-02-12
 */
public class SequenceAlignment implements Alignment {

	private static final DecimalFormat DFLONG = new DecimalFormat("0.0000");
	private static final DecimalFormat DFSHORT = new DecimalFormat("0.000");

	/**
	 * Sequence one; the top Sequence
	 */
	private final Sequence seq1;
	/**
	 * Sequence two; the bottom Sequence
	 */
	private final Sequence seq2;
	/**
	 * The two rows of the Alignment
	 */
	private final char[][] rows;
	/**
	 * The indices of the aligned residues, used by Paul format: map is an array
	 * of arrays whereas the first dimension just differs between sequence 1 or
	 * 2 the scd dimension is an array of length n or m referring to the length
	 * n or m of the sequences 1 or 2. every position in the scd dimension of
	 * the arrays contains the index of the residue from the other sequence it
	 * is aligned to
	 */
	private final int[][] map;
	/**
	 * The score of the Alignment
	 */
	private final double score;
	/**
	 * The length of the Alignment
	 */
	private final int length;

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
	public SequenceAlignment(Sequence seq1, Sequence seq2, char[][] rows,
			double score) {
		this.seq1 = seq1;
		this.seq2 = seq2;
		this.rows = rows;
		this.score = score;
		this.length = calcLength();
		map = calcMap();
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
	 * @param alignedResidues
	 *            the mapping from sequence to aligned sequence positions, see
	 *            javadoc of alignedResidues
	 */
	public SequenceAlignment(Sequence seq1, Sequence seq2, char[] row1,
			char[] row2, double score) {
		this.seq1 = seq1;
		this.seq2 = seq2;
		this.score = score;
		char[][] temp = new char[2][];
		temp[0] = row1;
		temp[1] = row2;
		this.rows = temp;
		this.length = row1.length;
		map = calcMap();
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
	 * @param alignedResidues
	 *            the mapping from sequence to aligned sequence positions, see
	 *            javadoc of alignedResidues
	 */
	public SequenceAlignment(Sequence seq1, Sequence seq2, String row1,
			String row2, double score) {
		this.seq1 = seq1;
		this.seq2 = seq2;
		this.score = score;
		char[][] temp = new char[2][];
		temp[0] = row1.toCharArray();
		temp[1] = row2.toCharArray();
		this.rows = temp;
		this.length = row1.length();
		map = calcMap();
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
	 * Calculates a 2D array similiar to the map[][], but that only contains
	 * the aligned residues
	 * 
	 * @return the 2D array
	 */
	public int[][] getAlignedResidues() {
		int[][] result = new int[2][countAlignedResidues()];
		int temp = 0;
		for (int i = 0; i < map[0].length; i++) {
			if (map[0][i] != -1) {
				result[0][temp] = i;
				result[1][temp] = map[0][i];
				temp++;
			}
		}
		return result;
	}

	/**
	 * Calculates the number of aligned residues
	 * 
	 * @return the number of aligned residues
	 */
	public int countAlignedResidues() {
		int result = 0;
		for (int i = 0; i < map[0].length; i++) {
			if (map[0][i] != -1) {
				result++;
			}
		}
		return result;
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
			if (rows[0][i] == '-') { // insertion
				result[1][col2] = -1;
				col2++;
			} else if (rows[1][i] == '-') { // deletion
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
	 * @return String of the format "id1 id2 score"
	 */
	public String toString() {
		String result = seq1.getID() + " " + seq2.getID() + " "
				+ DFLONG.format(score);
		return result;
	}

	/**
	 * 
	 * @return String of the format ">id1 id2 score\n id1: row1\n id2: row2"
	 */
	public String toStringVerbose() {
		String result = ">" + seq1.getID() + " " + seq2.getID() + " "
				+ DFSHORT.format(score) + "\n" + seq1.getID() + ": "
				+ String.valueOf(rows[0]) + "\n" + seq2.getID() + ": "
				+ String.valueOf(rows[1]);
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

	public Sequence getComponent(int n) {
		if (n == 0)
			return seq1;
		else
			return seq2;
	}

	public char[] getRow(int n) {
		return rows[n].clone();
	}

	public char[][] getRows() {
		return rows.clone();
	}

	public String getRowAsString(int n) {
		String res = "";
		for (char x : rows[n]) {
			res += x;
		}
		return res;
	}

	/**
	 * essentially copies the current alignment TODO this should be a
	 * constructor
	 */
	public SequenceAlignment duplicate() {
		return new SequenceAlignment(seq1, seq2, rows, score);
	}

	/**
	 * this neat function calculates the first / last aligned positions in an
	 * alignment (in both rows)
	 * 
	 * @author huberste
	 * @param ali
	 *            the alignment
	 * @return an int[4] with result[0] the number of the first aligned position
	 *         in the 1st sequence, result[1] the number of the last aligned
	 *         position in the 1st sequence, result[2] the number of the first
	 *         aligned position in the 2nd sequence and result[3] the number of
	 *         the last aligned position in the 2nd sequence.
	 */
	public static int[] calculateAlignedPositions(Alignment ali) {
		int[] result = new int[4];

		int[][] map = ali.calcMap();

		for (int i = 0; i < map[0].length; i++) {
			if (map[0][i] >= 0) {
				result[0] = i;
				break;
			}
		}

		for (int i = map[0].length - 1; i >= 0; i++) {
			if (map[0][i] >= 0) {
				result[1] = i;
				break;
			}
		}

		for (int i = 0; i < map[1].length; i++) {
			if (map[1][i] >= 0) {
				result[2] = i;
				break;
			}
		}

		for (int i = map[1].length - 1; i >= 0; i++) {
			if (map[1][i] >= 0) {
				result[3] = i;
				break;
			}
		}

		return result;
	}

}