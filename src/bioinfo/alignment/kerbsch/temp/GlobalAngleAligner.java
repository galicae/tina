package bioinfo.alignment.kerbsch.temp;

import java.io.BufferedWriter;
import java.io.IOException;
import java.io.OutputStreamWriter;

import bioinfo.Sequence;
import bioinfo.alignment.Alignable;
import bioinfo.alignment.Alignment;
import bioinfo.alignment.SequenceAlignment;
import bioinfo.alignment.gotoh.GlobalSequenceGotoh;
import bioinfo.alignment.gotoh.Gotoh;
import bioinfo.alignment.matrices.QuasarMatrix;
import bioinfo.proteins.PDBEntry;
import bioinfo.proteins.PDBFileReader;
import bioinfo.superpos.Kabsch;
import bioinfo.superpos.PDBReduce;
import bioinfo.superpos.TMMain;
import bioinfo.superpos.Transformation;

public class GlobalAngleAligner extends Gotoh {
	private static final int MM_ANGLEBORDER = (int) (8.9 * Gotoh.FACTOR);
	private static final int MP_ANGLEBORDER = (int) (16.9 * Gotoh.FACTOR);
	private static final int PM_ANGLEBORDER = (int) (99.26 * Gotoh.FACTOR);
	private static final int PP_ANGLEBORDER = (int) (58.69 * Gotoh.FACTOR);
	
	private static final int INIT_VAL = Integer.MIN_VALUE / 2;
	private int[][][][][] angleMatrix;
	private String dssppath = null;
	private int[] stances;

	public GlobalAngleAligner(double go, double ge, String anglepath,
			String dssppath) {
		super(go, ge);
		this.angleMatrix = DihedralAngles.readAngles(anglepath);
		this.dssppath = dssppath;
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
		stances = DihedralAngles.readStances(dssppath + "/"
				+ sequence1.getID() + ".dssp");

		for (int i = 2; i < sequence1.length(); i++) {
			for (int j = 2; j < sequence2.length(); j++) {
				tempScore[i - 1][j - 1] = getScore(seq1[i - 2]-65, seq1[i - 1]-65,
						seq1[i]-65, seq2[j - 2]-65, seq2[j - 1]-65, seq2[j]-65,
						stances[i - 1]);
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

		int x = seq1.length - 1;
		int y = seq2.length - 1;
		int score = M[x][y];
		String row0 = seq1[x] + "-";
		String row1 = "-" + seq2[y];
		int actScore = 0;
		char actx;
		char acty;
		while (x > 1 && y > 1) {

			actScore = M[x][y];
			actx = seq1[x - 1];
			acty = seq2[y - 1];

			if (actScore == M[x - 1][y - 1]
					+ getScore(seq1[x - 2]-65, actx-65, seq1[x]-65, seq2[y - 2]-65, acty-65,
							seq2[y]-65,stances[x-1])) {
				row0 += actx;
				row1 += acty;
				y--;
				x--;
			} else if (actScore == D[x][y]) {
				while (D[x][y] == D[x][y - 1] + gapExtend && y > 1) {
					row0 += "-";
					row1 += acty;
					y--;
					acty = seq2[y - 1];
				}
				row0 += "-";
				row1 += acty;
				y--;
			} else if (actScore == I[x][y]) {
				while (I[x][y] == I[x - 1][y] + gapExtend && x > 1) {
					row0 += actx;
					row1 += "-";
					x--;
					actx = seq1[x - 1];
				}
				row0 += actx;
				row1 += "-";
				x--;
			}
		}
		for (int i = y - 1; i >= 0; i--) {
			row0 += "-";
			row1 += seq2[i];
		}
		for (int i = x - 1; i >= 0; i--) {
			row0 += seq1[i];
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

	private int getScore(int pre1, int res1, int fol1, int pre2, int res2,
			int fol2, int stance) {
		int result;
		int border;
		
		if(stance == 0){
			border = MM_ANGLEBORDER;
		}else if (stance == 1){
			border = MP_ANGLEBORDER;
		}else if (stance == 2){
			border = PM_ANGLEBORDER;
		}else{
			border = PP_ANGLEBORDER;
		}
		result = border
				- ((Math.abs(angleMatrix[res1][pre1][fol1][0][stance]
						- angleMatrix[res2][pre2][fol2][0][stance]) + Math
							.abs(angleMatrix[res1][pre1][fol1][1][stance]
									- angleMatrix[res2][pre2][fol2][1][stance])) / 2);

		return result;
	}

//	public static void main(String[] args) throws Exception {
//		GlobalAngleAligner test = new GlobalAngleAligner(-20.0, -5.0, "angles","../GoBi_old/DSSP");
//		GlobalSequenceGotoh test2 = new GlobalSequenceGotoh(-12.0,-1.0,QuasarMatrix.DAYHOFF_MATRIX);
//
//		Sequence seq1 = new Sequence("1j2xA00","GPLDVQVTEDAVRRYLTRKPMTTKDLLKKFQTKKTGLSSEQTVNVLAQILKRLNPERKMINDKMHFSLK");
//		Sequence seq2 = new Sequence("1wq2B00","MEEAKQKVVDFLNSKSKSKFYFNDFTDLFPDMKQREVKKILTALVNDEVLEYWSSGSTTMYGLKG");
//		SequenceAlignment align1 = test.align(seq1, seq2);
//		SequenceAlignment align2 = test2.align(seq1, seq2);
//		System.out.println(align1.toStringVerbose());
//		System.out.println(align1.getScore()/ align1.countAlignedResidues());
//		
//		PDBFileReader pdbreader = new PDBFileReader("../GoBi_old/STRUCTURES");
//		TMMain superpos = new TMMain();
//		PDBEntry p = pdbreader.readFromFolderById("1j2xA00");
//		PDBEntry q = pdbreader.readFromFolderById("1wq2B00");
//		Transformation tr1 = superpos.calculateTransformation(align1, p, q);
//		Transformation tr2 = superpos.calculateTransformation(align2, p, q);
//		System.out.println("TMScore für Angles: " + tr1.getTmscore());
//		System.out.println("TMScore für Sequence: " + tr2.getTmscore());
//	}
}
