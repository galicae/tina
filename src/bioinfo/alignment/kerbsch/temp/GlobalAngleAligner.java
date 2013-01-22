package bioinfo.alignment.kerbsch.temp;

import java.util.HashMap;
import java.util.Map.Entry;

import bioinfo.Sequence;
import bioinfo.alignment.Alignable;
import bioinfo.alignment.Alignment;
import bioinfo.alignment.SequenceAlignment;
import bioinfo.alignment.gotoh.Gotoh;

public class GlobalAngleAligner extends Gotoh {
	private static final int INIT_VAL = Integer.MIN_VALUE / 2;
	private static final int ANGLE_BORDER = 3000;
	private HashMap<Character, HashMap<String, int[]>> angleprofiles = new HashMap<Character,HashMap<String,int[]>>();

	public GlobalAngleAligner(double go, double ge,
			HashMap<Character, HashMap<String, double[]>> angleprofiles) {
		super(go, ge);

		for (Entry<Character, HashMap<String, double[]>> amino : angleprofiles
				.entrySet()) {
			this.angleprofiles
					.put(amino.getKey(), new HashMap<String, int[]>());
			for (Entry<String, double[]> profile : amino.getValue().entrySet()) {
				this.angleprofiles.get(amino.getKey()).put(profile.getKey(),
						new int[2]);
				this.angleprofiles.get(amino.getKey()).get(profile.getKey())[0] = (int) (profile
						.getValue()[0] * Gotoh.FACTOR);
				this.angleprofiles.get(amino.getKey()).get(profile.getKey())[1] = (int) (profile
						.getValue()[1] * Gotoh.FACTOR);
			}
		}
	}

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

		for (int i = 2; i < sequence1.length(); i++) {
			for (int j = 2; j < sequence2.length(); j++) {
				tempScore[i - 1][j - 1] = score(seq1[i - 2], seq1[i - 1],
						seq1[i], seq2[j - 2], seq2[j - 1], seq2[j]);
			}
		}

		for (int i = 2; i < sequence1.length(); i++) {
			for (int j = 2; j < sequence2.length(); j++) {
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
		char[] seq1 = ((Sequence) sequence1).getSequence();
		char[] seq2 = ((Sequence) sequence2).getSequence();
		
		int x = sequence1.length() - 1;
		int y = sequence2.length() - 1;
		int score = M[x + 1][y + 1];
		String row0 = "";
		String row1 = "";
		int actScore = 0;
		char actx;
		char acty;
		while (x > 1 && y > 1) {

			actScore = M[x + 1][y + 1];
			actx = seq1[x-1];
			acty = seq2[y-1];

			if (actScore == M[x][y] + score(seq1[x-2],actx, seq1[x], seq2[y-2], acty, seq2[y])) {
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
	 * @param the
	 *            two indices of the actual considered aminos in the 2 Sequences
	 * @return the similarity between the angle-profiles score between two
	 *         components of Alignable
	 */
	private int score(char seq1_pre, char seq1_res, char seq1_fol,
			char seq2_pre, char seq2_res, char seq2_fol) {
		int phiDifference = Math.abs(this.angleprofiles.get(seq1_res).get(
				"" + seq1_pre + seq1_fol)[0]
				- this.angleprofiles.get(seq2_res)
						.get("" + seq2_pre + seq2_fol)[0]);
		int psiDifference = Math.abs(this.angleprofiles.get(seq1_res).get(
				"" + seq1_pre + seq1_fol)[1]
				- this.angleprofiles.get(seq2_res)
						.get("" + seq2_pre + seq2_fol)[1]);
		int phiscore = ANGLE_BORDER - phiDifference;
		int psiscore = ANGLE_BORDER - psiDifference;
		return phiscore+psiscore;
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
	
	public static void main(String[] args){
		GlobalAngleAligner test = new GlobalAngleAligner(-110.0,-102.0,DihedralAngles.read("angles"));
		
		Sequence seq1 = new Sequence("id1","MFKVYGYDSNIHKCVYCDNAKRLLTVKKQPFEFINIMPEKGVFDDEKIAELLTKLGRDTQIGLTMPQVFAPDGSHIGGFDQLREYF");
		Sequence seq2 = new Sequence("id2","MFKVYGYDSNIHKCVYCDNAKRLLTVKKQPFEFINIMPEKGVFDDEKIAELLTKLGRDTQIGLTMPQVFAPDGSHIGGFDQLREYF");
		SequenceAlignment alignment = test.align(seq1, seq2);
		System.out.println(alignment.toStringVerbose());
	}
}
