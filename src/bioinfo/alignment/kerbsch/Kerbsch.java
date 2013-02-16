package bioinfo.alignment.kerbsch;

import java.io.IOException;
import java.util.HashMap;

import bioinfo.Sequence;
import bioinfo.alignment.Alignable;
import bioinfo.alignment.Alignment;
import bioinfo.alignment.SequenceAlignment;
import bioinfo.alignment.gotoh.Gotoh;

public class Kerbsch extends Gotoh {

	private static final int INIT_VAL = Integer.MIN_VALUE / 2;

	// E-Value Params
	private final double EVALUE_CUTOFF_HB = 5.00;
	private final double EVALUE_CUTOFF_POL = 10.00;
	private final double EVALUE_CUTOFF_DAYHOFF = 6.00;
	private final double EVALUE_CUTOFF_SECSTRUCT = 9.0;

	// Lambda Params
	private final double LAMBDA_HB = 0.616;
	private final double LAMBDA_POL = 0.357;
	private final double LAMBDA_DAYHOFF = 0.38;
	private final double LAMBDA_SECSTRUCT = 0.408;

	private final int DBLENGTH = 63500;

	// size for matrices
	private int xsize;
	private int ysize;

	char[] sec1;
	char[] sec2;
	char[] seq1;
	char[] seq2;

	// gapopens, gapextends for substitution matrices
	private int hbGapOpen = -4;
	private int hbGapExtend = -2;

	private int polGapOpen = -10;
	private int polGapExtend = -1;

	private int secStructGapOpen = -4;
	private int secStructGapExtend = -3;

	private int seqGapOpen = -11;
	private int seqGapExtend = -1;

	// feature weights
	private final int hbWeight = (int) (0.1 * Gotoh.FACTOR);
	private final int polWeight = (int) (0.1 * Gotoh.FACTOR);
	private final int secStructWeight = (int) (0.3 * Gotoh.FACTOR);
	private final int seqWeight = (int) (0.1 * Gotoh.FACTOR);

	private FBGotoh hbLocals;
	private FBGotoh polLocals;
	private FBGotoh secStructLocals;
	private FBGotoh dayhoffLocals;
	private HashMap<String, char[]> seclib;

	// FBGotoh matrices
	private int[][] hbM;
	private int[][] polM;
	private int[][] secStructM;
	private int[][] dayhoffM;
	private int[][] tempscore;

	public Kerbsch(double gapOpen, double gapExtend, int[][] seqMatrix,
			int[][] hbMatrix, int[][] polMatrix, int[][] secStructMatrix,
			HashMap<String, char[]> seclib) throws IOException {
		super(gapOpen, gapExtend);

		hbLocals = new FBGotoh(hbGapOpen, hbGapExtend, EVALUE_CUTOFF_HB,
				LAMBDA_HB, DBLENGTH, hbMatrix);
		polLocals = new FBGotoh(polGapOpen, polGapExtend, EVALUE_CUTOFF_POL,
				LAMBDA_POL, DBLENGTH, polMatrix);
		secStructLocals = new FBGotoh(secStructGapOpen,
				EVALUE_CUTOFF_SECSTRUCT, LAMBDA_SECSTRUCT, DBLENGTH,
				secStructGapExtend, secStructMatrix);
		dayhoffLocals = new FBGotoh(seqGapOpen, seqGapExtend,
				EVALUE_CUTOFF_DAYHOFF, LAMBDA_DAYHOFF, DBLENGTH, seqMatrix);
		this.seclib = seclib;
	}

	public Kerbsch(double gapOpen, double gapExtend, double[][] seqMatrix,
			double[][] hbMatrix, double[][] polMatrix,
			double[][] secStructMatrix, HashMap<String, char[]> seclib) throws IOException {
		super(gapOpen, gapExtend);

		// for (int i = 0; i < seqMatrix.length; i++) {
		// for (int j = 0; j < seqMatrix[0].length; j++) {
		// seqMatrix[i][j] = (int) (Gotoh.FACTOR * seqMatrix[i][j]);
		// hbMatrix[i][j] = (int) (Gotoh.FACTOR * hbMatrix[i][j]);
		// polMatrix[i][j] = (int) (Gotoh.FACTOR * polMatrix[i][j]);
		// secStructMatrix[i][j] = (int) (Gotoh.FACTOR * secStructMatrix[i][j]);
		// }
		// }
		hbLocals = new FBGotoh(hbGapOpen, hbGapExtend, EVALUE_CUTOFF_HB,
				LAMBDA_HB, DBLENGTH, hbMatrix);
		polLocals = new FBGotoh(polGapOpen, polGapExtend, EVALUE_CUTOFF_POL,
				LAMBDA_POL, DBLENGTH, polMatrix);
		secStructLocals = new FBGotoh(secStructGapOpen,
				EVALUE_CUTOFF_SECSTRUCT, LAMBDA_SECSTRUCT, DBLENGTH,
				secStructGapExtend, secStructMatrix);
		dayhoffLocals = new FBGotoh(seqGapOpen, seqGapExtend,
				EVALUE_CUTOFF_DAYHOFF, LAMBDA_DAYHOFF, DBLENGTH, seqMatrix);
		this.seclib = seclib;
	}

	public SequenceAlignment align(Alignable sequence1, Alignable sequence2) {
		this.xsize = sequence1.length() + 1;
		this.ysize = sequence2.length() + 1;

		this.M = new int[xsize][ysize];
		this.I = new int[xsize][ysize];
		this.D = new int[xsize][ysize];
		tempscore = new int[xsize][ysize];
		
		this.sequence1 = (Sequence) sequence1;
		this.sequence2 = (Sequence) sequence2;
		this.seq1 = ((Sequence) sequence1).getSequence();
		this.seq2 = ((Sequence) sequence2).getSequence();
		this.sec1 = seclib.get(sequence1.getID());
		this.sec2 = seclib.get(sequence2.getID());

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
		hbM = hbLocals.align(seq1, seq2);
		polM = polLocals.align(seq1, seq2);
		secStructM = secStructLocals.align(sec1, sec2);
		dayhoffM = dayhoffLocals.align(seq1, seq2);

		// calc tempscore
		for (int i = 0; i < xsize; i++) {
			for (int j = 0; j < ysize; j++) {
				tempscore[i][j] = hbM[i][j] + polM[i][j] + secStructM[i][j]
						+ dayhoffM[i][j];
//				System.out.print(dayhoffM[i][j]+"\t");
			}
//			System.out.println();
		}

		for (int i = 1; i < xsize; i++) {
			for (int j = 1; j < ysize; j++) {

				// Kerbsch matrix
				D[i][j] = Math.max(M[i][j - 1] + gapOpen + gapExtend,
						D[i][j - 1] + gapExtend);
				I[i][j] = Math.max(M[i - 1][j] + gapOpen + gapExtend,
						I[i - 1][j] + gapExtend);
				M[i][j] = Math.max(M[i - 1][j - 1] + tempscore[i][j],
						Math.max(I[i][j], D[i][j]));
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
			row1 += seq2[i - 1];
		}
		for (int i = M.length - 1; i > x + 1; i--) {
			row0 += seq1[i - 1];
			row1 += "-";
		}

		while (x >= 0 && y >= 0) {

			actScore = M[x + 1][y + 1];
			actx = seq1[x];
			acty = seq2[y];

			if (actScore == M[x][y] + tempscore[x + 1][y + 1]) {
				row0 += actx;
				row1 += acty;
				y--;
				x--;
			} else if (actScore == D[x + 1][y + 1]) {
				while (D[x + 1][y + 1] == D[x + 1][y] + gapExtend && y > 0) {
					row0 += "-";
					row1 += acty;
					y--;
					acty = seq2[y];
				}
				row0 += "-";
				row1 += acty;
				y--;
			} else if (actScore == I[x + 1][y + 1]) {
				while (I[x + 1][y + 1] == I[x][y + 1] + gapExtend && x > 0) {
					row0 += actx;
					row1 += "-";
					x--;
					actx = seq1[x];
				}
				row0 += actx;
				row1 += "-";
				x--;
			}
		}
		for (int i = y + 1; i > 0; i--) {
			row0 += "-";
			row1 += seq2[i - 1];
		}
		for (int i = x + 1; i > 0; i--) {
			row0 += seq1[i - 1];
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
