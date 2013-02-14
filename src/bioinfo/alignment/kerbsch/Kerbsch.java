package bioinfo.alignment.kerbsch;

import java.util.HashMap;

import bioinfo.Sequence;
import bioinfo.alignment.Alignable;
import bioinfo.alignment.Alignment;
import bioinfo.alignment.SequenceAlignment;
import bioinfo.alignment.gotoh.Gotoh;

public class Kerbsch extends Gotoh {

	private static final int INIT_VAL = Integer.MIN_VALUE / 2;

	// size for matrices
	private int xsize;
	private int ysize;

	Sequence sec1;
	Sequence sec2;

	// gapopens, gapextends for substitution matrices
	private int hbGapOpen = -5 * Gotoh.FACTOR;
	private int hbGapExtend = -1 * Gotoh.FACTOR;

	private int polGapOpen = -11 * Gotoh.FACTOR;
	private int polGapExtend = -1 * Gotoh.FACTOR;

	private int secStructGapOpen = -15 * Gotoh.FACTOR;
	private int secStructGapExtend = -4 * Gotoh.FACTOR;

	private int seqGapOpen = -12 * Gotoh.FACTOR;
	private int seqGapExtend = -1 * Gotoh.FACTOR;

	// feature weights
	private final int hbWeight = (int) (0.1 * Gotoh.FACTOR);
	private final int polWeight = (int) (0.1 * Gotoh.FACTOR);
	private final int secStructWeight = (int) (0.3 * Gotoh.FACTOR);
	private final int seqWeight = (int) (0.1 * Gotoh.FACTOR);

	private FBGotoh hbCore;
	private FBGotoh polCore;
	private FBGotoh secStructCore;
	private FBGotoh seqCore;
	private HashMap<String, char[]> seclib;

	// gotoh matrices
	private int[][] hbM;
	private int[][] polM;
	private int[][] secStructM;
	private int[][] seqM;

	public Kerbsch(double gapOpen, double gapExtend, int[][] seqMatrix,
			int[][] hbMatrix, int[][] polMatrix, int[][] secStructMatrix,
			HashMap<String, char[]> seclib) {
		super(gapOpen, gapExtend);

//		hbCore = new FBGotoh(hbGapOpen, hbGapExtend, hbMatrix);
//		polCore = new FBGotoh(polGapOpen, polGapExtend, polMatrix);
//		secStructCore = new FBGotoh(secStructGapOpen, secStructGapExtend,
//				secStructMatrix);
//		seqCore = new FBGotoh(seqGapOpen, seqGapExtend, seqMatrix);
//		this.seclib = seclib;
	}

	public Kerbsch(double gapOpen, double gapExtend, double[][] seqMatrix,
			double[][] hbMatrix, double[][] polMatrix,
			double[][] secStructMatrix, HashMap<String, char[]> seclib) {
		super(gapOpen, gapExtend);

		for (int i = 0; i < seqMatrix.length; i++) {
			for (int j = 0; j < seqMatrix[0].length; j++) {
				seqMatrix[i][j] = (int) (Gotoh.FACTOR * seqMatrix[i][j]);
				hbMatrix[i][j] = (int) (Gotoh.FACTOR * hbMatrix[i][j]);
				polMatrix[i][j] = (int) (Gotoh.FACTOR * polMatrix[i][j]);
				secStructMatrix[i][j] = (int) (Gotoh.FACTOR * secStructMatrix[i][j]);
			}
		}

//		this.seclib = seclib;

//		hbCore = new FBGotoh(hbGapOpen, hbGapExtend, hbMatrix);
//		polCore = new FBGotoh(polGapOpen, polGapExtend, polMatrix);
//		secStructCore = new FBGotoh(secStructGapOpen, secStructGapExtend,
//				secStructMatrix);
//		seqCore = new FBGotoh(seqGapOpen, seqGapExtend, seqMatrix);
	}

	public SequenceAlignment align(Alignable sequence1, Alignable sequence2) {
		this.xsize = sequence1.length() + 1;
		this.ysize = sequence2.length() + 1;

		this.M = new int[xsize][ysize];
		this.I = new int[xsize][ysize];
		this.D = new int[xsize][ysize];

		this.sequence1 = (Sequence) sequence1;
		this.sequence2 = (Sequence) sequence2;
		this.sec1 = new Sequence(sequence1.getID(), seclib.get(sequence1
				.getID()));
		this.sec2 = new Sequence(sequence2.getID(), seclib.get(sequence2
				.getID()));

		prepareMatrices();
		calculateMatrices();
		return (SequenceAlignment) traceback();
	}

	private void prepareMatrices() {
		for (int i = 0; i < xsize; i++) {
			D[i][0] = INIT_VAL;
		}

		for (int i = 0; i < ysize; i++) {
			I[0][i] = INIT_VAL;
		}
	}

	private void calculateMatrices() {

		this.hbM = hbCore.align(sequence1, sequence2);
		this.polM = polCore.align(sequence1, sequence2);
		this.secStructM = secStructCore.align(sec1, sec2);
		this.seqM = seqCore.align(sequence1, sequence2);

		for (int i = 1; i < xsize; i++) {
			for (int j = 1; j < ysize; j++) {

				// Kerbsch matrix
				D[i][j] = Math.max(M[i][j - 1] + gapOpen + gapExtend,
						D[i][j - 1] + gapExtend);
				I[i][j] = Math.max(M[i - 1][j] + gapOpen + gapExtend,
						I[i - 1][j] + gapExtend);
				M[i][j] = Math.max(M[i - 1][j - 1] + (hbM[i][j] * hbWeight)
						+ (polM[i][j] * polWeight)
						+ (secStructM[i][j] * secStructWeight)
						+ (seqM[i][j] * seqWeight), Math.max(I[i][j], D[i][j]));
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

			if (actScore == M[x][y] + (hbM[x + 1][y + 1] * hbWeight) + (polM[x + 1][y + 1] * polWeight)
					+ (secStructM[x + 1][y + 1] * secStructWeight) + (seqM[x + 1][y + 1] * seqWeight)) {
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
