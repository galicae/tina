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
public class StructureAlignment implements Alignment{

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
	private final int score;
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
	public StructureAlignment(PDBEntry seq1, PDBEntry seq2, AminoAcid[][] rows, int score) {
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
	
	public StructureAlignment
		(PDBEntry seq1, PDBEntry seq2, AminoAcid[] row1, AminoAcid[] row2, int score)
	{
		this.seq1 = seq1;
		this.seq2 = seq2;
		this.score = score;
		AminoAcid[][] temp = new AminoAcid[2][];
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
	public StructureAlignment
		(PDBEntry seq1, PDBEntry seq2, List<AminoAcid> row1, List<AminoAcid> row2, int score)
	{
		this.seq1 = seq1;
		this.seq2 = seq2;
		this.score=score;
		AminoAcid[][] temp = new AminoAcid[2][];
		temp[0] = row1.toArray(new AminoAcid[row1.size()]);
		temp[1] = row2.toArray(new AminoAcid[row2.size()]);
		this.rows = temp;
		this.length = row1.size();
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
			if (rows[0][i] == null) {	// insertion
				result[1][col2]=-1;
				col2++;
			} else if (rows[1][i] == null) {	// deletion
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
		String amino1 = "";
		String amino2 = "";
		
		for(int i = 0; i != rows[0].length; i++){
			if(rows[0][i] == null){
				amino1 += "-";
			}else{
				amino1 += rows[0][i].toString();
			}
			if(rows[1][i] == null){
				amino2 += "-";
			}else{
				amino2 += rows[1][i].toString();
			}
		}
		
		String result =
			">"+seq1.getID()+" "+seq2.getID()+" "+DFSHORT.format(score)+"\n"+
				seq1.getID() + ": " + amino1 + "\n"+
				seq2.getID() + ": " + amino2;
		return result;
	}
	
	/**
	 * 
	 * @return length of the Alignment
	 */
	public int length() {
		return length;
	}
	
	public int getScore() {
		return score;
	}
	
	public PDBEntry getComponent(int n) {
		if (n==0) return seq1;
		else return seq2;
	}
	
	public AminoAcid[] getRow(int n) {
		return rows[n];
	}
}
