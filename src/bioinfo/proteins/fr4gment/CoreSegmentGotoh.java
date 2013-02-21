package bioinfo.proteins.fr4gment;

import java.util.LinkedList;

import bioinfo.Sequence;
import bioinfo.alignment.Alignable;
import bioinfo.alignment.Alignment;
import bioinfo.alignment.SequenceAlignment;
import bioinfo.alignment.gotoh.Gotoh;
import bioinfo.proteins.fragm3nt.ProteinFragment;

/**
 * Local Alignment of two Sequences
 * 
 * @author andreseitz, galicae
 */
public class CoreSegmentGotoh extends Gotoh {

	private static final int INIT_VAL = Integer.MIN_VALUE / 2;
	private int[][] scoringmatrix;
	private ProteinFragment used;
	private ProteinFragment xFrag;
	private double cutoff;

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
	public CoreSegmentGotoh(double gapOpen, double gapExtend, double cutoff,
			ProteinFragment used, ProteinFragment x) {
		super(gapOpen, gapExtend);
		this.cutoff = cutoff;
		this.used = used;
		this.xFrag = x;
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
		score = 0;
		int begin = -1;
		int end = -1;

		int i = 0;
		while (begin == -1) {
			if (row0[i] != '-' && row1[i] != '-') {
				begin = i;
			}
			i++;
		}

		i = row0.length - 1;
		while (end == -1) {
			if (row0[i] != '-' && row1[i] != '-') {
				end = i;
			}
			i--;
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
	 * prepares matrices for local/global alignment
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
	 * calculates matrices using given scoring function and gap penalty
	 * 
	 */
	private void calculateMatrices() {
		scoringmatrix = new int[sequence1.length()][sequence2.length()];

		for (int i = 1; i <= sequence1.length(); i++) {
			for (int j = 1; j <= sequence2.length(); j++) {
				scoringmatrix[i - 1][j - 1] = score(i - 1, j - 1);
			}
		}
		enhance(scoringmatrix);

		for (int i = 1; i <= sequence1.length(); i++) {
			for (int j = 1; j <= sequence2.length(); j++) {
				D[i][j] = Math.max(M[i][j - 1] + gapOpen + gapExtend,
						D[i][j - 1] + gapExtend);
				I[i][j] = Math.max(M[i - 1][j] + 10 * gapOpen + 10 * gapExtend,
						I[i - 1][j] + 10 * gapExtend);
				M[i][j] = Math.max(M[i - 1][j - 1] + scoringmatrix[i - 1][j - 1],
						Math.max(I[i][j], Math.max(D[i][j], 0)));

			}
		}
	}

	private void enhance(int[][] tempScore) {
		LinkedList<int[]> posPoints = new LinkedList<int[]>();
		for (int i = 0; i < tempScore.length; i++) {
			for (int j = 0; j < tempScore[0].length; j++) {
				if (tempScore[i][j] >= 0) {
					int[] positive = { i, j };
					posPoints.add(positive);
				}
			}
		}

		for (int[] cur : posPoints) {
			int x = cur[0];
			int y = cur[1];
			if (x > 0) {
				tempScore[x - 1][y] = tempScore[x][y];
				if (y > 0)
					tempScore[x - 1][y - 1] = tempScore[x][y];
				if (y < tempScore[0].length - 1)
					tempScore[x - 1][y + 1] = tempScore[x][y];
			}
			if (x < tempScore.length - 1) {
				tempScore[x + 1][y] = tempScore[x][y];
				if (y > 0)
					tempScore[x + 1][y - 1] = tempScore[x][y];
				if (y < tempScore[0].length - 1)
					tempScore[x + 1][y + 1] = tempScore[x][y];
			}
		}
	}

	private Alignment traceback() {
		return null;
	}

	/**
	 * this traceback starts from a given column in the alignment in order to
	 * identify the
	 * 
	 * @return Alignment of the two given Alignables
	 */
	public int[] traceback(int start, int end) {

//		int start = pos[0];
//		int end = pos[1];
		int max = INIT_VAL;
		int x = 0;
		int y = 0;
		int resStart = 0;
		int resEnd = 0;
		for (int i = end; i < end + 1; i++) {
			for (int j = 0; j != M[i].length; j++) {
				if (max <= M[i][j]) {
					max = M[i][j];
					x = i - 1;
					y = j - 1;
				}
			}
		}
		resStart = y + 1;
		int score = max;
		String row0 = "";
		String row1 = "";
		int actScore = 0;
		char actx;
		char acty;

		while (x >= start && y >= 0 && M[x + 1][y + 1] != 0) {
//			System.out.println(x + " " + y);
			actScore = M[x + 1][y + 1];
			actx = (Character) sequence1.getComp(x);
			acty = (Character) sequence2.getComp(y);
			if (actScore == M[x][y] + scoringmatrix[x][y]) {
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
				while (I[x + 1][y + 1] == I[x][y + 1] + 10 * gapExtend && x > 0) {
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

		int[] result = { y, resStart };
		return result;
	}

	/**
	 * @param two
	 *            components of Alignable implementing equals
	 * @return score between two components of Alignable
	 */
	private int score(int x, int y) {
		double result = euclideanDistance(used.getResidue(x),
				this.xFrag.getResidue(y));
		result = 10* cutoff - result;
		return (int) result * Gotoh.FACTOR;
	}

	private double euclideanDistance(double[] x, double[] y) {
		double result = 0;
		for (int i = 0; i < x.length; i++) {
			result += Math.pow(x[i] - y[i], 2);
		}
		return Math.sqrt(result);
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
