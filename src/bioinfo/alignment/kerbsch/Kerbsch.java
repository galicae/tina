package bioinfo.alignment.kerbsch;

import bioinfo.Sequence;
import bioinfo.alignment.Alignable;
import bioinfo.alignment.Alignment;
import bioinfo.alignment.SequenceAlignment;
import bioinfo.alignment.gotoh.Gotoh;

public class Kerbsch extends Gotoh {

	private static final int INIT_VAL = Integer.MIN_VALUE / 2;

	// substitution matrices
	private int[][] hbMatrix;
	private int[][] polMatrix;
	private int[][] secStructMatrix;
	private int[][] substMatrix;

	// gapopens, gapextends for substitution matrices
	private int hbGapOpen = -400;
	private int polGapOpen = -400;
	private int secStructGapOpen = -100;
	private int substGapOpen = -1500;
	private int hbGapExtend = -100;
	private int polGapExtend = -100;
	private int secStructGapExtend = -50;
	private int substGapExtend = -300;

	// gotoh matrices
	private int[][] hbM;
	private int[][] hbD;
	private int[][] hbI;

	private int[][] polM;
	private int[][] polD;
	private int[][] polI;

	private int[][] secStructM;
	private int[][] secStructD;
	private int[][] secStructI;

	private int[][] substM;
	private int[][] substD;
	private int[][] substI;

	public Kerbsch(double gapOpen, double gapExtend, int[][] substMatrix,
			int[][] hbMatrix, int[][] polMatrix, int[][] secStructMatrix) {
		super(gapOpen, gapExtend);
		this.hbMatrix = hbMatrix;
		this.polMatrix = polMatrix;
		this.secStructMatrix = secStructMatrix;
		this.substMatrix = substMatrix;
	}

	public Kerbsch(double gapOpen, double gapExtend, double[][] substMatrix,
			double[][] hbMatrix, double[][] polMatrix,
			double[][] secStructMatrix) {
		super(gapOpen, gapExtend);
		this.hbMatrix = new int[hbMatrix.length][hbMatrix[0].length];
		this.polMatrix = new int[polMatrix.length][polMatrix[0].length];
		this.secStructMatrix = new int[secStructMatrix.length][secStructMatrix[0].length];
		this.substMatrix = new int[substMatrix.length][substMatrix[0].length];
		for (int i = 0; i < substMatrix.length; i++) {
			for (int j = 0; j < substMatrix[0].length; j++) {
				this.substMatrix[i][j] = (int) (Gotoh.FACTOR * substMatrix[i][j]);
				this.hbMatrix[i][j] = (int) (Gotoh.FACTOR * substMatrix[i][j]);
				this.polMatrix[i][j] = (int) (Gotoh.FACTOR * substMatrix[i][j]);
				this.secStructMatrix[i][j] = (int) (Gotoh.FACTOR * substMatrix[i][j]);
			}
		}
	}

	public SequenceAlignment align(Alignable sequence1, Alignable sequence2) {
		this.M = new int[sequence1.length() + 1][sequence2.length() + 1];
		this.I = new int[sequence1.length() + 1][sequence2.length() + 1];
		this.D = new int[sequence1.length() + 1][sequence2.length() + 1];

		this.hbM = new int[sequence1.length() + 1][sequence2.length() + 1];
		this.hbI = new int[sequence1.length() + 1][sequence2.length() + 1];
		this.hbD = new int[sequence1.length() + 1][sequence2.length() + 1];

		this.polM = new int[sequence1.length() + 1][sequence2.length() + 1];
		this.polI = new int[sequence1.length() + 1][sequence2.length() + 1];
		this.polD = new int[sequence1.length() + 1][sequence2.length() + 1];

		this.secStructM = new int[sequence1.length() + 1][sequence2.length() + 1];
		this.secStructI = new int[sequence1.length() + 1][sequence2.length() + 1];
		this.secStructD = new int[sequence1.length() + 1][sequence2.length() + 1];

		this.substM = new int[sequence1.length() + 1][sequence2.length() + 1];
		this.substI = new int[sequence1.length() + 1][sequence2.length() + 1];
		this.substD = new int[sequence1.length() + 1][sequence2.length() + 1];

		this.sequence1 = (Sequence) sequence1;
		this.sequence2 = (Sequence) sequence2;
		prepareMatrices();
		calculateMatrices();
		return (SequenceAlignment) traceback();
	}

	private void prepareMatrices() {
		for (int i = 0; i < sequence1.length(); i++) {
			D[i][0] = INIT_VAL;
			hbD[i][0] = INIT_VAL;
			polD[i][0] = INIT_VAL;
			secStructD[i][0] = INIT_VAL;
			substD[i][0] = INIT_VAL;
		}

		for (int i = 0; i < sequence2.length(); i++) {
			I[0][i] = INIT_VAL;
			hbI[0][i] = INIT_VAL;
			polI[0][i] = INIT_VAL;
			secStructI[0][i] = INIT_VAL;
			substI[0][i] = INIT_VAL;
		}
	}

	private int score(int[][] matrix, char x, char y) {
		return matrix[x - 65][y - 65];
	}

	private void calculateMatrices() {
		char[] seq1 = ((Sequence) sequence1).getSequence();
		char[] seq2 = ((Sequence) sequence2).getSequence();

		for (int i = 1; i <= sequence1.length(); i++) {
			for (int j = 1; j <= sequence2.length(); j++) {
				// substitution gotoh matrices
				hbD[i][j] = Math.max(hbM[i][j - 1] + hbGapOpen + hbGapExtend,
						hbD[i][j - 1] + hbGapExtend);
				hbI[i][j] = Math.max(hbM[i - 1][j] + hbGapOpen + hbGapExtend,
						hbI[i - 1][j] + hbGapExtend);
				hbM[i][j] = Math.max(
						hbM[i - 1][j - 1]
								+ score(hbMatrix, seq1[i - 1], seq2[j - 1]),
						Math.max(hbI[i][j], hbD[i][j]));
				polD[i][j] = Math.max(polM[i][j - 1] + polGapOpen
						+ polGapExtend, polD[i][j - 1] + polGapExtend);
				polI[i][j] = Math.max(polM[i - 1][j] + polGapOpen
						+ polGapExtend, polI[i - 1][j] + polGapExtend);
				polM[i][j] = Math.max(
						polM[i - 1][j - 1]
								+ score(polMatrix, seq1[i - 1], seq2[j - 1]),
						Math.max(polI[i][j], polD[i][j]));
				secStructD[i][j] = Math.max(secStructM[i][j - 1]
						+ secStructGapOpen + secStructGapExtend,
						secStructD[i][j - 1] + secStructGapExtend);
				secStructI[i][j] = Math.max(secStructM[i - 1][j]
						+ secStructGapOpen + secStructGapExtend,
						secStructI[i - 1][j] + secStructGapExtend);
				secStructM[i][j] = Math.max(
						secStructM[i - 1][j - 1]
								+ score(secStructMatrix, seq1[i - 1],
										seq2[j - 1]),
						Math.max(secStructI[i][j], secStructD[i][j]));
				substD[i][j] = Math.max(substM[i][j - 1] + substGapOpen
						+ substGapExtend, substD[i][j - 1] + substGapExtend);
				substI[i][j] = Math.max(substM[i - 1][j] + substGapOpen
						+ substGapExtend, substI[i - 1][j] + substGapExtend);
				substM[i][j] = Math.max(
						substM[i - 1][j - 1]
								+ score(substMatrix, seq1[i - 1], seq2[j - 1]),
						Math.max(substI[i][j], substD[i][j]));

				// Kerbsch matrix
				D[i][j] = Math.max(M[i][j - 1] + gapOpen + gapExtend,
						D[i][j - 1] + gapExtend);
				I[i][j] = Math.max(M[i - 1][j] + gapOpen + gapExtend,
						I[i - 1][j] + gapExtend);
				M[i][j] = Math.max(M[i - 1][j - 1] + hbMatrix[i][j]
						+ polMatrix[i][j] + secStructMatrix[i][j]
						+ substMatrix[i][j], Math.max(I[i][j], D[i][j]));
			}
		}
	}

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

			if (actScore == M[x][y] + hbMatrix[x][y] + polMatrix[x][y]
					+ secStructMatrix[x][y] + substMatrix[x][y]) {
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
	 * flips a char[] on itself
	 * 
	 * @param in
	 *            the character array in question
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
