package bioinfo.alignment.kerbsch;

import java.util.HashMap;

import bioinfo.Sequence;
import bioinfo.alignment.Alignable;
import bioinfo.alignment.Alignment;
import bioinfo.alignment.SequenceAlignment;
import bioinfo.alignment.gotoh.Gotoh;

public class FreeshiftMusterLite extends Gotoh {

	private static final int INIT_VAL = Integer.MIN_VALUE / 2;

	// size for matrices
	private int xsize;
	private int ysize;

	// substitution matrices
	private int[][] hbMatrix = new int[26][26];
	private int[][] polMatrix = new int[26][26];
	private int[][] secStructMatrix = new int[26][26];
	private int[][] substMatrix = new int[26][26];
	
	//feature weights
	private final int hbWeight;
	private final int polWeight;
	private final int secStructWeight;
	private final int substWeight;

	private int[][] tempScore;
	private HashMap<String,char[]> secStructSeqlib;
	
	public FreeshiftMusterLite(double gapOpen, double gapExtend, HashMap<String,char[]> sslib, int[][] hbMatrix,
			int[][] polMatrix, int[][] secStructMatrix, int[][] substMatrix, int hbWeight, int polWeight, int ssWeight, int substWeight) {
		super(gapOpen, gapExtend);
		
		this.gapOpen *= Gotoh.FACTOR;
		this.gapExtend *= Gotoh.FACTOR;
		
		this.hbMatrix = hbMatrix;
		this.polMatrix = polMatrix;
		this.secStructMatrix = secStructMatrix;
		this.substMatrix = substMatrix;
		this.hbWeight = hbWeight;
		this.polWeight = polWeight;
		this.secStructWeight = ssWeight;
		this.substWeight = substWeight;
		this.secStructSeqlib = sslib;
	}

	public FreeshiftMusterLite(double gapOpen, double gapExtend, HashMap<String,char[]> sslib, double[][] hbMatrix,
			double[][] polMatrix, double[][] secStructMatrix,
			double[][] substMatrix, double hbWeight, double polWeight, double ssWeight ,double substWeight) {
		super(gapOpen, gapExtend);
				
		for (int i = 0; i < substMatrix.length; i++) {
			for (int j = 0; j < substMatrix[i].length; j++) {
				this.substMatrix[i][j] = (int) (Gotoh.FACTOR * substMatrix[i][j]);
				this.hbMatrix[i][j] = (int) (Gotoh.FACTOR * hbMatrix[i][j]);
				this.polMatrix[i][j] = (int) (Gotoh.FACTOR * polMatrix[i][j]);
				this.secStructMatrix[i][j] = (int) (Gotoh.FACTOR * secStructMatrix[i][j]);
			}
		}
		
		this.hbWeight = (int)(hbWeight * Gotoh.FACTOR);
		this.polWeight = (int)(polWeight * Gotoh.FACTOR);
		this.secStructWeight = (int)(ssWeight * Gotoh.FACTOR);
		this.substWeight = (int)(substWeight * Gotoh.FACTOR);
		this.secStructSeqlib = sslib;
	}

	public SequenceAlignment align(Alignable sequence1, Alignable sequence2) {
		this.xsize = sequence1.length() + 1;
		this.ysize = sequence2.length() + 1;

		this.M = new int[xsize][ysize];
		this.I = new int[xsize][ysize];
		this.D = new int[xsize][ysize];
		tempScore = new int[xsize][ysize];

		this.sequence1 = (Sequence) sequence1;
		this.sequence2 = (Sequence) sequence2;
		prepareMatrices();
		calculateMatrices();
		return (SequenceAlignment) traceback();
	}

	private void prepareMatrices() {
		for (int i = 1; i < xsize; i++) {
			D[i][0] = INIT_VAL;
		}

		for (int i = 1; i < ysize; i++) {
			I[0][i] = INIT_VAL;
		}
	}

	private int score(int[][] matrix, char x, char y) {
		return matrix[x - 65][y - 65];
	}

	private void calculateMatrices() {	
		char[] seq1 = ((Sequence) sequence1).getSequence();
		char[] seq2 = ((Sequence) sequence2).getSequence();
		
		char[] sec1 = secStructSeqlib.get(sequence1.getId());
		char[] sec2 = secStructSeqlib.get(sequence2.getId());
		
		for (int i = 1; i < xsize; i++) {
			for (int j = 1; j < ysize; j++) {
				tempScore[i][j] = score(hbMatrix, seq1[i - 1], seq2[j - 1]) * hbWeight
						+ score(polMatrix, seq1[i - 1], seq2[j - 1]) * polWeight
						+ score(secStructMatrix, sec1[i - 1], sec2[j - 1]) * secStructWeight
						+ score(substMatrix, seq1[i - 1], seq2[j - 1]) * substWeight;
			}
		}	

		for (int i = 1; i < xsize; i++) {
			for (int j = 1; j < ysize; j++) {
				D[i][j] = Math.max(M[i][j - 1] + gapOpen + gapExtend,
						D[i][j - 1] + gapExtend);
				I[i][j] = Math.max(M[i - 1][j] + gapOpen + gapExtend,
						I[i - 1][j] + gapExtend);
				M[i][j] = Math.max(M[i - 1][j - 1] + tempScore[i][j],
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

			if (actScore == M[x][y] + tempScore[x+1][y+1]) {
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
