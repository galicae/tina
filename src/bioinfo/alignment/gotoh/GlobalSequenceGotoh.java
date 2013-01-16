package bioinfo.alignment.gotoh;

import java.util.ArrayList;
import java.util.List;

import bioinfo.Sequence;
import bioinfo.alignment.Alignable;
import bioinfo.alignment.Alignment;
import bioinfo.alignment.SequenceAlignment;

/**
 * Global Alignment of two sequences
 * 
 * @author andreseitz
 */
public class GlobalSequenceGotoh extends Gotoh {

	private static final int INIT_VAL = Integer.MIN_VALUE / 2;
	private int[][] scoringmatrix;

	// Sequence sequence1;
	// Sequence sequence2;

	/**
	 * 
	 * @param gapOpen
	 * @param gapExtend
	 * @param scoringmatrix
	 *            26x26 matrix containing all scoring values plus some empty
	 *            lines for faster access
	 */
	public GlobalSequenceGotoh(double gapOpen, double gapExtend,
			int[][] scoringmatrix) {
		super(gapOpen, gapExtend);
		this.scoringmatrix = scoringmatrix;
	}

	public GlobalSequenceGotoh(double gapOpen, double gapExtend,
			double[][] scoringmatrix) {
		super(gapOpen, gapExtend);
		this.scoringmatrix = new int[scoringmatrix.length][scoringmatrix[0].length];
		for (int i = 0; i != scoringmatrix.length; i++) {
			for (int j = 0; j != scoringmatrix[0].length; j++) {
				this.scoringmatrix[i][j] = (int) (Gotoh.FACTOR * scoringmatrix[i][j]);
			}
		}
	}

	@Override
	public SequenceAlignment align(Alignable sequence1, Alignable sequence2) {
		this.M = new int[sequence1.length() + 1][sequence2.length() + 1];
		this.I = new int[sequence1.length() + 1][sequence2.length() + 1];
		this.D = new int[sequence1.length() + 1][sequence2.length() + 1];
		this.sequence1 = sequence1;
		this.sequence2 = sequence2;
		prepareMatrices();
		calculateMatrices();
		return (SequenceAlignment) traceback();
	}

	@Override
	public boolean check(Alignment alignment) {
		SequenceAlignment ali = (SequenceAlignment) alignment;
		int score = 0;
		char[] row0 = ali.getRow(0);
		char[] row1 = ali.getRow(1);
		for (int i = 0; i < row0.length; i++) {
			if (row0[i] == '-' || row1[i] == '-') {
				score += gapOpen;
				while (i < row0.length && (row0[i] == '-' || row1[i] == '-')) {
					score += this.gapExtend;
					i++;
				}
				i--;
			} else {
				score += score(row0[i], row1[i]);
			}
		}
		if (1.0d * score / Gotoh.FACTOR == ali.getScore()) {
			return true;
		} else {
			return false;
		}
	}

	/**
	 * prepares matrices for global alignment
	 * 
	 */
	private void prepareMatrices() {
		for (int i = 1; i <= sequence1.length(); i++) { // old: hSeq
			// I[i][0] = 0; //old: vGap
			D[i][0] = INIT_VAL; // old: hGap
			if (i == 1) {
				M[i][0] = M[i - 1][0] + gapOpen + gapExtend;
			} else {
				M[i][0] = M[i - 1][0] + gapExtend;
			}
		}

		for (int i = 1; i <= sequence2.length(); i++) { // old: vSeq
			// D[0][i] = 0;
			if (i == 1) {
				M[0][i] = M[0][i - 1] + gapOpen + gapExtend;
			} else {
				M[0][i] = M[0][i - 1] + gapExtend;
			}
			I[0][i] = INIT_VAL;
		}
		D[0][0] = INIT_VAL;
		I[0][0] = INIT_VAL;
		M[0][0] = 0;
	}

	/**
	 * calculates matrices using given scoring function and gap penalty
	 * 
	 */
	private void calculateMatrices() {
		int[][] tempScore = new int[sequence1.length()][sequence2.length()];
		char[] seq1 = ((Sequence) sequence1).getSequence();
		char[] seq2 = ((Sequence) sequence2).getSequence();

		for (int i = 1; i <= sequence1.length(); i++) {
			for (int j = 1; j <= sequence2.length(); j++) {
				tempScore[i - 1][j - 1] = score(seq1[i - 1], seq2[j - 1]);
			}
		}

		for (int i = 1; i <= sequence1.length(); i++) {
			for (int j = 1; j <= sequence2.length(); j++) {
				D[i][j] = Math.max(M[i][j - 1] + gapOpen + gapExtend,
						D[i][j - 1] + gapExtend);
				I[i][j] = Math.max(M[i - 1][j] + gapOpen + gapExtend,
						I[i - 1][j] + gapExtend);
				M[i][j] = Math.max(M[i - 1][j - 1] + tempScore[i - 1][j - 1],
						Math.max(I[i][j], D[i][j]));

			}
		}
	}

	/**
	 * Override this method in extensions!
	 * 
	 * @return Alignment of the two given Alignables
	 */
	private Alignment traceback() {
		List<int[]> map = new ArrayList<int[]>();
		
		int x = sequence1.length() - 1;
		int y = sequence2.length() - 1;
		int score = M[x + 1][y + 1];
		String row0 = "";
		String row1 = "";
		int actScore = 0;
		char actx;
		char acty;
		while (x >= 0 && y >= 0) {

			actScore = M[x + 1][y + 1];
			actx = (Character) sequence1.getComp(x);
			acty = (Character) sequence2.getComp(y);

			if (actScore == M[x][y] + score(actx, acty)) {
				row0 += actx;
				row1 += acty;
				map.add(new int[]{x,y}); //store aligned indices of the two sequences
				y--;
				x--;
			} else if (actScore == D[x + 1][y + 1]) {
				while (D[x + 1][y + 1] == D[x + 1][y] + gapExtend && y > 0) {
					row0 += "-";
					row1 += acty;
					y--;
					acty = (Character) sequence2.getComp(y);
				}
				row0 += "-";
				row1 += acty;
				y--;
			} else if (actScore == I[x + 1][y + 1]) {
				while (I[x + 1][y + 1] == I[x][y + 1] + gapExtend && x > 0) {
					row0 += actx;
					row1 += "-";
					x--;
					actx = (Character) sequence1.getComp(x);
				}
				row0 += actx;
				row1 += "-";
				x--;
			}
		}
		for (int i = y + 1; i > 0; i--) {
			row0 += "-";
			row1 += sequence2.getComp(i - 1);
		}
		for (int i = x + 1; i > 0; i--) {
			row0 += sequence1.getComp(i - 1);
			row1 += "-";
		}

		return new SequenceAlignment((Sequence) sequence1,
				(Sequence) sequence2, flip(row0.toCharArray()),
				flip(row1.toCharArray()), 1.0d * score / Gotoh.FACTOR, map.toArray(new int[map.size()][]));
	}

	/**
	 * @param two
	 *            components of Alignable implementing equals
	 * @return score between two components of Alignable
	 */
	private int score(char x, char y) {
		return scoringmatrix[x - 65][y - 65];
	}

	/**
	 * flips a char[] on itself
	 * @param in the character array in question
	 * @return the reversed array
	 */
	private char[] flip(char[] in) {
		char[] out = new char[in.length];
		for (int i = in.length - 1; i >= 0; i--) {
			out[out.length - 1 - i] = in[i];
		}
		return out;
	}

}
