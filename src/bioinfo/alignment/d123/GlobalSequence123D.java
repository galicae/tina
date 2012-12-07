package bioinfo.alignment.d123;

import bioinfo.Sequence;
import bioinfo.alignment.Alignable;
import bioinfo.alignment.Alignment;
import bioinfo.alignment.SequenceAlignment;
import bioinfo.alignment.gotoh.Gotoh;
import bioinfo.proteins.SSCCEntry;

/**
 * 
 * @author gruppe_4 Freeshift version of 123D
 */
public class GlobalSequence123D extends Gotoh {

	private static final int INIT_VAL = Integer.MIN_VALUE / 2;
	private int[] secStruct, localConts, globalConts;
	private int[][] scoringmatrix, secStrucPref, weights;
	private int[][][] contactPot;

	/**
	 * 
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
	 * @param secondaryStructurePreferences
	 *            a matrix containing secondary structure preference for all
	 *            amino acids
	 * @param contactPot
	 *            the contact potentials of the amino acids
	 * @param sscc
	 *            an sscc entry containing the structural information concerning
	 *            the template structure
	 */
	public GlobalSequence123D(double[][] scoringmatrix,
			double[][] secondaryStructurePreferences, double[][] weights,
			double[][][] contactPot, SSCCEntry sscc) {
		super(0.0d, 0.0d);

		this.secStruct = sscc.getSecStruct();
		this.localConts = sscc.getLocalConts();
		this.globalConts = sscc.getGlobalConts();

		this.contactPot = new int[contactPot.length][contactPot[0].length][contactPot[0][0].length];
		for (int i = 0; i != contactPot.length; i++) {
			for (int j = 0; j != contactPot[0].length; j++) {
				for (int k = 0; k != contactPot[0][0].length; k++) {
					this.contactPot[i][j][k] = (int) (Gotoh.FACTOR * contactPot[i][j][k]);
				}
			}
		}

		this.secStrucPref = new int[secondaryStructurePreferences.length][secondaryStructurePreferences[0].length];
		for (int i = 0; i != secondaryStructurePreferences.length; i++) {
			for (int j = 0; j != secondaryStructurePreferences[0].length; j++) {
				this.secStrucPref[i][j] = (int) (Gotoh.FACTOR * secondaryStructurePreferences[i][j]);
			}
		}
		this.weights = new int[weights.length][weights[0].length];
		for (int i = 0; i != weights.length; i++) {
			for (int j = 0; j != weights[0].length; j++) {
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
	 * function that calculates the score of a match between amino acids x and y
	 * given the secondary structure of y
	 * 
	 * @param x
	 *            amino acid x
	 * @param y
	 *            amino acid y
	 * @param stY
	 *            the secondary structure of y
	 * @return the score of x matching y
	 */
	private int match(char x, char y, int stY) {
		int seqScore = score(x, y);
		int prefScore = secStrucPref[stY][x - 1];
		int lcontScore = contactPot[stY][x - 1][localConts[y - 1]];
		int gcontScore = contactPot[stY][x - 1][globalConts[y - 1]];
		int result = weights[4][stY] * lcontScore + weights[5][stY]
				* gcontScore + weights[3][stY] * prefScore + weights[1][stY]
				* seqScore;
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
		char[] seq1 = ((Sequence) sequence1).getSequence();
		char[] seq2 = ((Sequence) sequence2).getSequence();

		int[] gapOpen = new int[3];
		gapOpen[0] = this.gapOpen * weights[1][0];
		gapOpen[1] = this.gapOpen * weights[1][1];
		gapOpen[2] = this.gapOpen * weights[1][2];

		int[] gapExtend = new int[3];
		gapExtend[0] = this.gapExtend * weights[2][0];
		gapExtend[1] = this.gapExtend * weights[2][1];
		gapExtend[2] = this.gapExtend * weights[2][2];

		int[][][] tempScore = new int[sequence1.length()][sequence2.length()][3];
		for (int i = 1; i <= sequence1.length(); i++) {
			for (int j = 1; j <= sequence2.length(); j++) {
				int strY = secStruct[j - 1];
				tempScore[i][j][strY] = match(seq1[i - 1], seq2[j - 1], strY);
			}
		}

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
				score += gapOpen[secStruct[row1[i]]];
				while (i <= end && (row0[i] == '-' || row1[i] == '-')) {
					score += gapExtend[secStruct[row1[i]]];
					i++;
				}
				i--;
			} else {
				score += tempScore[row0[i]][row1[i]][secStruct[row1[0]]];
			}
		}
		if (1.0d * score / (Gotoh.FACTOR * Gotoh.FACTOR)== ali.getScore()) {
			return true;
		} else {
			return false;
		}
	}

	/**
	 * prepares matrices for freeshift alignment
	 * 
	 */
	private void prepareMatrices() {
		int[] gapOpen = new int[3];
		gapOpen[0] = this.gapOpen * weights[1][0];
		gapOpen[1] = this.gapOpen * weights[1][1];
		gapOpen[2] = this.gapOpen * weights[1][2];

		int[] gapExtend = new int[3];
		gapExtend[0] = this.gapExtend * weights[2][0];
		gapExtend[1] = this.gapExtend * weights[2][1];
		gapExtend[2] = this.gapExtend * weights[2][2];
		
		for (int i = 1; i <= sequence1.length(); i++) { // old: hSeq
			// I[i][0] = 0; //old: vGap
			D[i][0] = INIT_VAL; // old: hGap
			if (i == 1) {
				M[i][0] = M[i - 1][0] + gapOpen[secStruct[0]] + gapExtend[secStruct[0]];
			} else {
				M[i][0] = M[i - 1][0] + gapExtend[secStruct[0]];
			}
		}

		for (int i = 1; i <= sequence2.length(); i++) { // old: vSeq
			// D[0][i] = 0;
			if (i == 1) {
				M[0][i] = M[0][i - 1] + gapOpen[secStruct[i-1]] + gapExtend[secStruct[i-1]];
			} else {
				M[0][i] = M[0][i - 1] + gapExtend[secStruct[i-1]];
			}
			I[0][i] = INIT_VAL;
		}
		D[0][0] = INIT_VAL;
		I[0][0] = INIT_VAL;
		M[0][0] = 0;
	}

	/**
	 * calculates matrices using scoring-function and gap-penalty
	 * 
	 */
	private void calculateMatrices() {
		// getting everything out of the loop as in Gotoh; it seemed to help A
		// LOT
		char[] seq1 = ((Sequence) sequence1).getSequence();
		char[] seq2 = ((Sequence) sequence2).getSequence();

		int[] gapOpen = new int[3];
		gapOpen[0] = this.gapOpen * weights[1][0];
		gapOpen[1] = this.gapOpen * weights[1][1];
		gapOpen[2] = this.gapOpen * weights[1][2];

		int[] gapExtend = new int[3];
		gapExtend[0] = this.gapExtend * weights[2][0];
		gapExtend[1] = this.gapExtend * weights[2][1];
		gapExtend[2] = this.gapExtend * weights[2][2];

		int[][][] tempScore = new int[sequence1.length()][sequence2.length()][3];
		for (int i = 1; i <= sequence1.length(); i++) {
			for (int j = 1; j <= sequence2.length(); j++) {
				int strY = secStruct[j - 1];
				tempScore[i][j][strY] = match(seq1[i - 1], seq2[j - 1], strY);
			}
		}

		// now the main loop where stuff is actually computed
		for (int i = 1; i <= sequence1.length(); i++) {
			for (int j = 1; j <= sequence2.length(); j++) {
				int strY = secStruct[j - 1];
				D[i][j] = Math.max(M[i][j - 1] + gapOpen[strY]
						+ gapExtend[strY], D[i][j - 1] + gapExtend[strY]);
				I[i][j] = Math.max(M[i - 1][j] + gapOpen[strY]
						+ gapExtend[strY], I[i - 1][j] + gapExtend[strY]);
				M[i][j] = Math.max(M[i - 1][j - 1]
						+ tempScore[i - 1][j - 1][strY],
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

		int x = sequence1.length() - 1;
		int y = sequence2.length() - 1;
		int score = M[x + 1][y + 1];
		String row0 = "";
		String row1 = "";
		int actScore = 0;
		int strY = -1;
		char actx;
		char acty;
		
		int[] gapOpen = new int[3];
		gapOpen[0] = this.gapOpen * weights[1][0];
		gapOpen[1] = this.gapOpen * weights[1][1];
		gapOpen[2] = this.gapOpen * weights[1][2];

		int[] gapExtend = new int[3];
		gapExtend[0] = this.gapExtend * weights[2][0];
		gapExtend[1] = this.gapExtend * weights[2][1];
		gapExtend[2] = this.gapExtend * weights[2][2];
		
		while (x >= 0 && y >= 0) {
			
			actScore = M[x + 1][y + 1];
			actx = (Character) sequence1.getComp(x);
			acty = (Character) sequence2.getComp(y);
			strY = secStruct[acty];

			if (actScore == M[x][y] + score(actx, acty)) {
				row0 += actx;
				row1 += acty;
				y--;
				x--;
			} else if (actScore == D[x + 1][y + 1]) {
				while (D[x + 1][y + 1] == D[x + 1][y] + gapExtend[strY] && y > 0) {
					row0 += "-";
					row1 += acty;
					y--;
					acty = (Character) sequence2.getComp(y);
				}
				row0 += "-";
				row1 += acty;
				y--;
			} else if (actScore == I[x + 1][y + 1]) {
				while (I[x + 1][y + 1] == I[x][y + 1] + gapExtend[strY] && x > 0) {
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
