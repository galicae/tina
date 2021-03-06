/******************************************************************************
 * bioinfo.alignment.Threading.java                                           *
 *                                                                            *
 * Contains the Class Threading which is an implementation of a               *
 * Structure - Sequence Alignment.                                            *
 *                                                                            *
 * This file is best read at line width 80 and tab width 4.                   *
 *                                                                   huberste *
 ******************************************************************************/
package bioinfo.alignment;

import java.text.DecimalFormat;

import bioinfo.Sequence;
import bioinfo.proteins.PDBEntry;

/**
 * A Threading is an Alignment of a (template) Structure and a (target)
 * Sequence.
 * 
 * @author huberste
 * @lastchange 2013-02-20
 */
public class Threading implements Alignment {

	/**
	 * template
	 */
	private PDBEntry structure;
	/**
	 * target
	 */
	private Sequence sequence;

	/**
	 * 2D array with the positions of row[0] template and row[1] target. <br />
	 * row[0][n] = -1 means that the n-th position in the alignment is an
	 * insertion and row[1][n] means that the n-th position in the alignment is
	 * an deletion. e.g.: <br />
	 * 012-3456-789 <br />
	 * -012-34-567-
	 */
	private int[][] rows;

	/**
	 * Score of the Threading.
	 */
	private double score;

	/**
	 * 2D array with map[0] the template and map[1] the target. <br />
	 * map[x][n] = -1 means that the n-th residue of the template (x=0) or
	 * target (x=1) are aligned to a gap and every other value means they are
	 * aligned to that value of the other Alignable.
	 */
	private int[][] map;

	/**
	 * positionsInAlignables[i][pos] is the position of the (pos-th character in
	 * the rows) in the i-th Alignable
	 */
	private int[][] positionsInAlignables;

	/**
	 * The length of each row of the Threading.
	 */
	private int length;

	/**
	 * Standard constructor for a new Threading
	 * 
	 * @param structure
	 * @param sequence
	 * @param rows
	 * @param score
	 */
	public Threading(PDBEntry structure, Sequence sequence, int[][] rows,
			double score) {
		this.structure = structure;
		this.sequence = sequence;
		this.setRows(rows);
		this.score = score;
	}

	/**
	 * sets the rows and updates the map, length and other stuff.
	 * 
	 * @param rows
	 * @throws Exception
	 */
	private void setRows(int[][] rows) {
		// Error Handling
		// if (rows.length != 2) {
		// throw new Exception("rows must be of length 2!");
		// }
		// if (rows[0].length != rows[1].length) {
		// throw new Exception("rows[0] and rows[1] must be the same length!");
		// }
		// set rows
		if (rows == null) {
			rows = new int[2][structure.length() + sequence.length()];
			for (int i = 0; i < structure.length(); i++) {
				rows[0][i] = i;
				rows[1][i] = -1;
			}
			for (int i = 0; i < sequence.length(); i++) {
				rows[0][structure.length() + i] = -1;
				rows[1][structure.length() + i] = i;
			}
		}
		this.rows = rows;
		// update length
		this.length = rows[0].length;
		// update Map
		this.map = calcMap();
		// update positionsInAlignables
		this.positionsInAlignables = calcPositions();

	}

	/**
	 * Calculates a 2D array similiar to the map[][], but that only contains the
	 * aligned residues
	 * 
	 * @return the 2D array
	 */
	@Override
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
	 * calculates the map between seq1 and seq2.<br />
	 * map[i][j] is -1 if j-th character of row(i) is aligned to a gap and the
	 * number of the character of row(1-i) it is aligned to else.
	 * 
	 * @return the map between seq1 and seq2
	 */
	@Override
	public int[][] calcMap() {
		int[][] result = new int[2][];
		result[0] = new int[structure.length()];
		result[1] = new int[sequence.length()];
		int col1 = 0;
		int col2 = 0;
		for (int i = 0; i < rows[0].length; i++) {
			if (rows[0][i] == -1) { // insertion
				result[1][col2] = -1;
				col2++;
			} else if (rows[1][i] == -1) { // deletion
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
	 * @return
	 */
	public int[][] calcPositions() {
		int[][] result = new int[2][rows[0].length];
		int structpos = 0;
		int seqpos = 0;
		for (int i = 0; i < rows[0].length; i++) {
			if (rows[0][i] == -1) { // insertion
				result[0][i] = -1;
				result[1][i] = seqpos;
				seqpos++;
			} else if (rows[1][i] == -1) { // deletion
				result[0][i] = structpos;
				result[1][i] = -1;
				structpos++;
			} else { // match
				result[0][i] = structpos;
				result[1][i] = seqpos;
				structpos++;
				seqpos++;
			}

		}
		return result;
	}

	/**
	 * @return the length of the Threading
	 */
	@Override
	public int length() {
		return length;
	}

	/**
	 * 
	 * @param score
	 */
	public void setScore(double score) {
		this.score = score;
	}

	/**
	 * @return the score of the Threading
	 */
	@Override
	public double getScore() {
		return score;
	}

	public int[][] getRows() {
		return rows;
	}

	/**
	 * 
	 * @return the rows as a String representation
	 */
	public String[] getRowsAsString() {
		String[] result = new String[2];
		result[0] = "";
		result[1] = "";
		int temppos = 0;
		int targpos = 0;
		for (int i = 0; i < rows[0].length; i++) {
			if (rows[0][i] != -1) {
				result[0] += structure.getComp(temppos).getName()
						.getOneLetterCode();
				temppos++;
			} else
				result[0] += '-';
			if (rows[1][i] != -1) {
				result[1] += sequence.getComp(targpos).toString();
				targpos++;
			} else
				result[1] += '-';
		}

		return result;
	}

	/**
	 * 
	 * @return the rows as a char Array representation
	 */
	public char[][] getRowsAsCharArray() {
		char[][] result = new char[2][rows[0].length];
		int temppos = 0;
		int targpos = 0;
		for (int i = 0; i < rows[0].length; i++) {
			if (rows[0][i] != -1) {
				result[0][i] = structure.getComp(temppos).getName()
						.getOneLetterCode().charAt(0);
				temppos++;
			} else
				result[0][i] = '-';
			if (rows[1][i] != -1) {
				result[1][i] = sequence.getComp(targpos);
				targpos++;
			} else
				result[1][i] = '-';
		}

		return result;
	}

	@Override
	public Alignable getComponent(int number) {
		return number == 0 ? structure : sequence;
	}

	@Override
	public Alignment duplicate() {
		Alignment result = null;
		try {
			result = new Threading(this.structure, this.sequence, this.rows,
					this.score);
		} catch (Exception e) {
			// never happens :-)
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
	 * 
	 * @param n
	 * @return the position od the n-th Residue in the Rows in the Structure
	 */
	public int getPositionInStructure(int n) {
		return positionsInAlignables[0][n];
	}

	public int getFirstAfterInStructure(int n) {
		while (n < positionsInAlignables[0].length
				&& positionsInAlignables[0][n] == -1) {
			n++;
		}
		return positionsInAlignables[0][n];
	}

	/**
	 * 
	 * @param n
	 * @return the position od the n-th Residue in the Rows in the Sequence
	 */
	public int getPositionInSequence(int n) {
		return positionsInAlignables[1][n];
	}

	public int getFirstAfterInSequence(int n) {
		while (n < positionsInAlignables[1].length
				&& positionsInAlignables[1][n] == -1) {
			n++;
		}
		return positionsInAlignables[1][n];
	}

	public Sequence getSequence() {
		return sequence;
	}

	public PDBEntry getStructure() {
		return structure;
	}

	public SequenceAlignment asSequenceAlignment() {
		SequenceAlignment result = new SequenceAlignment(
				structure.getSequence(), sequence, getRowsAsCharArray(), score);
		return result;
	}

	@Override
	public String toString() {
		return toStringVerbose();
		// DecimalFormat df = new DecimalFormat("0.0000");
		// return "Threading: structure: " + structure.getLongID()
		// + ", sequence: " + sequence.getID() + ", score: "
		// + df.format(score) + "\n";
	}

	@Override
	public String toStringVerbose() {
		DecimalFormat df = new DecimalFormat("0.0000");
		String result = "> " + "Threading: structure: " + structure.getID()
				+ ", sequence: " + sequence.getId() + ", score: "
				+ df.format(score) + "\n";
		String[] temp = getRowsAsString();
		result += structure.getID() + ": " + temp[0] + "\n";
		result += sequence.getId() + ": " + temp[1] + "\n";
		return result;
	}

}

/******************************************************************************
 * "A question that sometimes drives me hazy:                                 *
 *  Am I or are the others crazy?"                                            *
 * - Albert Einstein (1879 - 1955)                                            *
 ******************************************************************************/
