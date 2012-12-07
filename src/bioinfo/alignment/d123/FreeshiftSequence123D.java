package bioinfo.alignment.d123;

import bioinfo.Sequence;
import bioinfo.alignment.Alignable;
import bioinfo.alignment.Alignment;
import bioinfo.alignment.SequenceAlignment;
import bioinfo.alignment.gotoh.Gotoh;

/**
 * 
 * @author gruppe_4 Freeshift version of 123D
 */
public class FreeshiftSequence123D extends Gotoh {

	private static final int INIT_VAL = Integer.MIN_VALUE / 2;
	private int[][] scoringmatrix;
	private int[][] weights;
	private int[][] secStrucPref;
	private int[] secStruct;

	/**
	 * @param weights
	 *            the weights for sequence, secondary structure and contact
	 *            potential scores
	 * @param gapOpen
	 *            the sequence gap open penalty
	 * @param gapExtend
	 *            the sequence gap extend penalty
	 * @param scoringmatrix
	 *            26x26 matrix containing all scoring values plus some empty
	 *            lines for faster access
	 */
	public FreeshiftSequence123D(double[][] scoringmatrix,
			double[][] secondaryStructurePreferences, double[][] weights) {
		super(0.0d, 0.0d);
		this.secStrucPref = new int[secondaryStructurePreferences.length][secondaryStructurePreferences[0].length];
		for(int i = 0; i != secondaryStructurePreferences.length; i++) {
			for(int j = 0; j != secondaryStructurePreferences[0].length; j++) {
				this.secStrucPref[i][j] = (int) (Gotoh.FACTOR * secondaryStructurePreferences[i][j]);
			}
		}
		this.weights = new int[weights.length][weights[0].length];
		for(int i = 0; i != weights.length; i++) {
			for(int j = 0; j != weights[0].length; j++) {
				this.weights[i][j] = (int) (Gotoh.FACTOR * weights[i][j]);
			}
		}
		this.scoringmatrix = new int[scoringmatrix.length][scoringmatrix[0].length];
		for (int i = 0; i != scoringmatrix.length; i++) {
			for (int j = 0; j != scoringmatrix[0].length; j++) {
				this.scoringmatrix[i][j] = (int) (Gotoh.FACTOR * scoringmatrix[i][j]);
			}
		}
	}
	
	
	/**
	 * function that calculates the score of a match between amino acids x and y given the secondary structure of y
	 * @param x amino acid x
	 * @param y amino acid y
	 * @param stY the secondary structure of y
	 * @return the score of x matching y
	 */
	private double match(char x, char y, int stY) {
		//TODO: find some way to incorporate contacts information
		int seqScore = score(x, y);
		int prefScore = secStrucPref[stY][x - 1];
		int lcontScore = contactPot[stY][x - 1][localConts[y - 1]];
		int gcontScore = contactPot[stY][x - 1][globalConts[y - 1]];
		int result = weights[4][stY] * lcontScore + 
		weights[5][stY] * gcontScore + 
		weights[3][stY] * prefScore + 
		weights[1][stY] * seqScore;
		return result;
	}

	@Override
	public SequenceAlignment align(Alignable sequence1, Alignable sequence2) {
		this.M = new int[sequence1.length() + 1][sequence2.length() + 1];
		this.I = new int[sequence1.length() + 1][sequence2.length() + 1];
		this.D = new int[sequence1.length() + 1][sequence2.length() + 1];
		this.sequence1 = (Sequence) sequence1;
		this.sequence2 = (Sequence) sequence2;
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
		int begin = -1;
		int end = -1;

		int i = 0;
		if (row0[i] == '-') {
			while (row0[i] == '-') {
				i++;
			}
			begin = i;
		} else if (row1[0] == '-') {
			while (row1[i] == '-') {
				i++;
			}
			begin = i;
		} else {
			begin = i;
		}

		i = row1.length - 1;
		if (row1[i] == '-') {
			while (row1[i] == '-') {
				i--;
			}
			end = i;
		} else if (row0[i] == '-') {
			while (row0[i] == '-') {
				i--;
			}
			end = i;
		} else {
			end = i;
		}

		for (i = begin; i <= end; i++) {
			if (row0[i] == '-' || row1[i] == '-') {
				score += this.gapOpen;
				while (i <= end && (row0[i] == '-' || row1[i] == '-')) {
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
			// M[i][0] = 0;
		}

		for (int i = 1; i <= sequence2.length(); i++) { // old: vSeq
			// D[0][i] = 0;
			// M[0][i] = 0;
			I[0][i] = INIT_VAL;
		}
		D[0][0] = INIT_VAL;
		I[0][0] = INIT_VAL;
		// M[0][0] = 0;
	}

	/**
	 * calculates matrices using scoring-function and gap-penalty
	 * 
	 */
	private void calculateMatrices() {
		for (int i = 1; i <= sequence1.length(); i++) {
			for (int j = 1; j <= sequence2.length(); j++) {
				D[i][j] = Math.max(M[i][j - 1] + gapOpen + gapExtend,
						D[i][j - 1] + gapExtend);
				// System.out.println((M[i][j-1]+gapOpen+gapExtend)+" "+(D[i][j-1]+gapExtend));
				I[i][j] = Math.max(M[i - 1][j] + gapOpen + gapExtend,
						I[i - 1][j] + gapExtend);
				M[i][j] = Math.max(
						M[i - 1][j - 1]
								+ score((Character) sequence1.getComp(i - 1),
										(Character) sequence2.getComp(j - 1)),
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

		int max = INIT_VAL;
		int x = 0;
		int y = 0;

		for (int i = 0; i != M.length; i++) {
			if (M[i][M[i].length - 1] >= max) {
				max = M[i][M[i].length - 1];
				x = i - 1;
				y = M[i].length - 2;
			}
		}
		for (int i = (M[M.length - 1].length - 1); i >= 0; i--) {
			if (M[(M.length - 1)][i] > max) {
				max = M[(M.length - 1)][i];
				y = i - 1;
				x = M.length - 2;
			}
		}

		int score = max;
		String row0 = "";
		String row1 = "";
		int actScore = 0;
		char actx;
		char acty;

		for (int i = M[M.length - 1].length - 1; i > y + 1; i--) {
			row0 += "-";
			row1 += sequence2.getComp(i - 1);
		}
		for (int i = M.length - 1; i > x + 1; i--) {
			row0 += sequence1.getComp(i - 1);
			row1 += "-";
		}

		while (x >= 0 && y >= 0) {

			actScore = M[x + 1][y + 1];
			actx = (Character) sequence1.getComp(x);
			acty = (Character) sequence2.getComp(y);

			if (actScore == M[x][y] + score(actx, acty)) {
				row0 += actx;
				row1 += acty;
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
				flip(row1.toCharArray()), 1.0d * score / Gotoh.FACTOR);
	}

	/**
	 * @param two
	 *            components of Alignable implementing equals
	 * @return score between two components of Alignable
	 */
	private int score(char x, char y) {
		return scoringmatrix[x - 65][y - 65];
	}

	private char[] flip(char[] in) {
		char[] out = new char[in.length];
		for (int i = in.length - 1; i >= 0; i--) {
			out[out.length - 1 - i] = in[i];
		}
		return out;
	}

}