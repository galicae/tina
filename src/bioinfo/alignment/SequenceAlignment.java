/******************************************************************************
 * bioinfo.alignment.Alignment                                                *
 *                                                                            *
 * This file is best read at line width 80 and tab width 4.                   *
 ******************************************************************************/
package bioinfo.alignment;

import java.text.DecimalFormat;
import java.util.Locale;

import bioinfo.Sequence;

/**
 * @author gobi_4
 * @date November 25, 2012
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
	 * The score of the Alignment
	 */
	private final double score;
	/**
	 * The length of the Alignment
	 */
	private final int length;
	
	/**
	 * Constructs an Alignment of seq1 and seq2 with the given map and score.
	 * @param seq1 the 1st Sequence
	 * @param seq2 the 2nd Sequence
	 * @param rows the two rows of the alignment
	 * @param score the score of the Alignment
	 */
	public SequenceAlignment(Sequence seq1, Sequence seq2, char[][] rows, double score) {
		this.seq1 = seq1;
		this.seq2 = seq2;
		this.rows = rows;
		this.score=score;
		this.length = calcLength();
	}
	
	/**
	 * Constructs an Alignment of seq1 and seq2 with the given alignment rows
	 * and score.
	 * @param seq1 the 1st Sequence
	 * @param seq2 the 2nd Sequence
	 * @param row1 the first row of the Alignment
	 * @param row2 the second row of the Alignment
	 * @param score the score of the Alignment
	 */
	public SequenceAlignment
		(Sequence seq1, Sequence seq2, char[] row1, char[] row2, double score)
	{	
		this.seq1 = seq1;
		this.seq2 = seq2;
		this.score = score;
		char[][] temp = new char[2][];
		temp[0] = row1;
		temp[1] = row2;
		this.rows = temp;
		this.length = row1.length;
	}
	
	/**
	 * Constructs an Alignment of seq1 and seq2 with the given alignment rows
	 * and score.
	 * @param seq1 the 1st Sequence
	 * @param seq2 the 2nd Sequence
	 * @param row1 the first row of the Alignment
	 * @param row2 the second row of the Alignment
	 * @param score the score of the Alignment
	 */
	public SequenceAlignment
		(Sequence seq1, Sequence seq2, String row1, String row2, double score)
	{	
		this.seq1 = seq1;
		this.seq2 = seq2;
		this.score=score;
		char[][] temp = new char[2][];
		temp[0] = row1.toCharArray();
		temp[1] = row2.toCharArray();
		this.rows = temp;
		this.length = row1.length();
	}
	
	/**
	 * Calculates the length of the given alignment.
	 * @return the length of the alignment
	 */
	private int calcLength() {
		return rows[0].length;
	}
	
	/**
	 * calculates the map between seq1 and seq2.<br />
	 * map[i][j] is -1 if j-th character of row(i) is aligned to a gap and the
	 * number of the character of row(1-i) it is aligned to else.
	 * @return the map between seq1 and seq2
	 */
	public int[][] calcMap() {
		int[][] result = new int[2][];
		result[0] = new int[seq1.length()];
		result[1] = new int[seq2.length()];
		int col1 = 0;
		int col2 = 0;
		for (int i = 0; i < seq1.length(); i++) {
			if (rows[0][i] == '-') {	// insertion
				result[1][col2]=-1;
				col2++;
			} else if (rows[1][i] == '-') {	// deletion
				result[0][col1]=-1;
				col1++;
			} else {	// (mis)match
				result[0][col1]=col2;
				result[1][col2]=col1;
				col1++; col2++;
			}
		}
		return result;
	}
	
	/**
	 * 
	 * @return String of the format "idone idtwo score"
	 */
	public String toString() {
		String result =
				seq1.getID() + " " + seq2.getID() + " " + DFLONG.format(score);
		return result;
	}
	
	/**
	 * 
	 * @return String of the format ">idone idtwo score\n idone: rowone\n
	 *  idtwo: rowtwo"
	 */
	public String toStringVerbose() {
		String result =
			">"+seq1.getID()+" "+seq2.getID()+" "+DFSHORT.format(score)+"\n"+
				seq1.getID() + ": " + String.valueOf(rows[0]) + "\n"+
				seq2.getID() + ": " + String.valueOf(rows[1]);
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
		if (n==0) return seq1;
		else return seq2;
	}
	
	public char[] getRow(int n) {
		return rows[n];
	}
	
	public String getRowAsString(int n){
		String res = "";
		for(char x: rows[n]){
			res += x;
		}
		return res;
	}
	
}