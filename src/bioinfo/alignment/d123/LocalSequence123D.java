package bioinfo.alignment.d123;

import bioinfo.Sequence;
import bioinfo.alignment.Alignable;
import bioinfo.alignment.Alignment;
import bioinfo.alignment.SequenceAlignment;
import bioinfo.alignment.gotoh.Gotoh;
import bioinfo.proteins.SSCCEntry;
import bioinfo.proteins.SSCCLine;

/**
 * 
 * @author gruppe_4 Local version of 123D
 */
public class LocalSequence123D extends Gotoh {

	private static final int INIT_VAL = Integer.MIN_VALUE / 2;
	private int[] secStruct, localConts, globalConts, gapOpen, gapExtend = null;
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
	public LocalSequence123D(double go, double ge, double[][] scoringmatrix,
			double[][] secondaryStructurePreferences, double[][] weights,
			double[][][] contactPot) {
		super(go, ge);
		
		go = go*Gotoh.FACTOR;
		ge = ge*Gotoh.FACTOR;
		
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
		gapOpen = new int[3];
		gapOpen[0] = (int)(go * this.weights[1][0]);
		gapOpen[1] = (int)(go * this.weights[1][1]);
		gapOpen[2] = (int)(go * this.weights[1][2]);

		gapExtend = new int[3];
		gapExtend[0] = (int)(ge * this.weights[2][0]);
		gapExtend[1] = (int)(ge * this.weights[2][1]);
		gapExtend[2] = (int)(ge * this.weights[2][2]);
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
//		System.out.println("x: "+x);
//		System.out.println("y: "+y);
		int seqScore = score(x, y);
		int prefScore = secStrucPref[stY][x - 65];
		int lcontScore = contactPot[stY][localConts[y - 65]][x - 65];
		int gcontScore = contactPot[stY][globalConts[y - 65]][x - 65];
		int result = (weights[4][stY] * lcontScore) + 
					(weights[5][stY] * gcontScore) + 
					(weights[3][stY] * prefScore) + 
					(weights[0][stY] * seqScore);
//		System.out.println(stY+"      "+weights[0][stY]+" "+seqScore+" "+(seqScore*weights[0][stY])+"        "+weights[3][stY]+" "+prefScore+" "+(prefScore*weights[3][stY])+"       "+weights[4][stY]+" "+lcontScore+" "+(weights[4][stY]*lcontScore)+"       "+weights[5][stY]+" "+gcontScore+" "+(weights[5][stY]*gcontScore));
//		System.out.println(seqScore + " "+ weights[1][stY]);
//		System.out.println(prefScore+ " "+ weights[3][stY]);
//		System.out.println(lcontScore+ " "+ weights[4][stY]);
//		System.out.println(gcontScore+ " "+ weights[5][stY]);
//		System.out.println(result);
		return result;
	}

	/**
	 * this function precedes the actual alignment; here the SSCC file is read
	 * and the secondary structure saved in secStruct.
	 * 
	 * @param sequence1
	 *            the query sequence
	 * @param sequence2
	 *            the template sequence
	 * @param sscc
	 *            the SSCC file, containing secondary structure and contacts
	 *            information concerning the template structure
	 * @return the result of the align() function
	 */
	public SequenceAlignment align(Alignable sequence1, Alignable sequence2, SSCCEntry sscc){
		//parse SSCCEntry-----------------------------
		this.globalConts = new int[sscc.length()];
		this.localConts = new int[sscc.length()];
		this.secStruct = new int[sscc.length()];
		int localtemp;
		int globaltemp;
		
		SSCCLine sscctemp;
		for (int i = 0; i < sscc.length(); i++) {
			sscctemp = (SSCCLine)sscc.getComp(i);
			
			switch (sscctemp.getSecStruct()){
				case 'a': this.secStruct[i] = 0;break;
				case 'b': this.secStruct[i] = 1;break;
				case 'o': this.secStruct[i] = 2;break;
				default: this.secStruct[i] = 2;
			}
			localtemp = sscctemp.getLocCont();
			globaltemp = sscctemp.getGlobCont();
			
			if(localtemp > 13){
				this.localConts[i] = 13;
			}else{
				this.localConts[i] = sscctemp.getLocCont();
			}
			if(globaltemp > 13){
				this.globalConts[i] = 13;
			}else{
				this.globalConts[i] = sscctemp.getGlobCont();
			}	
		}
		//----------------------------------------------
		return align(sequence1, sequence2);
	}
	
	
	//caution: if calling align function without SSCCEntry then return is null
	@Override
	public SequenceAlignment align(Alignable sequence1, Alignable sequence2) {
		if(this.localConts == null || this.globalConts == null || this.secStruct == null){
			return null;
		}
		else{
			this.M = new int[sequence1.length() + 1][sequence2.length() + 1];
			this.I = new int[sequence1.length() + 1][sequence2.length() + 1];
			this.D = new int[sequence1.length() + 1][sequence2.length() + 1];
			this.sequence1 = (Sequence) sequence1;
			this.sequence2 = (Sequence) sequence2;
			prepareMatrices();
			calculateMatrices();
			return (SequenceAlignment) traceback();
		}
	}

	@Override
	public boolean check(Alignment alignment) {
		char[] seq1 = ((Sequence) sequence1).getSequence();
		char[] seq2 = ((Sequence) sequence2).getSequence();

		int[][] tempScore = new int[sequence1.length()][sequence2.length()];
		for (int i = 1; i <= sequence1.length(); i++) {
			for (int j = 1; j <= sequence2.length(); j++) {
				int strY = secStruct[j - 1];
				tempScore[i][j] = match(seq1[i - 1], seq2[j - 1], strY);
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
				score += tempScore[row0[i]][row1[i]];
			}
		}
		if (1.0d * score / (Gotoh.FACTOR * Gotoh.FACTOR)== ali.getScore()) {
			return true;
		} else {
			return false;
		}
	}

	/**
	 * prepares matrices for local alignment
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
	 * calculates matrices using the given scoring function and gap penalty
	 * 
	 */
	private void calculateMatrices() {

		// getting everything out of the loop as in Gotoh; it seemed to help A
		// LOT
		char[] seq1 = ((Sequence) sequence1).getSequence();
		char[] seq2 = ((Sequence) sequence2).getSequence();

		int[][] tempScore = new int[sequence1.length()][sequence2.length()];
		int strY;
		for (int i = 0; i < sequence1.length(); i++) {
			//System.out.println();
			for (int j = 0; j < sequence2.length(); j++) {
				strY = secStruct[j];
				tempScore[i][j] = match(seq1[i], seq2[j], strY);
				//System.out.print(String.format("%8.3f",(tempScore[i][j]/1000000.0d))+"\t");
			}
		}

		// now the main loop where stuff is actually computed
		for (int i = 1; i <= sequence1.length(); i++) {
			for (int j = 1; j <= sequence2.length(); j++) {
				strY = secStruct[j-1];
				D[i][j] = Math.max(M[i][j-1] + gapOpen[strY] + gapExtend[strY], D[i][j-1] + gapExtend[strY]);
				I[i][j] = Math.max(M[i-1][j] + gapOpen[strY] + gapExtend[strY], I[i-1][j] + gapExtend[strY]);
				M[i][j] = Math.max(M[i-1][j-1] + tempScore[i-1][j-1], Math.max(I[i][j], Math.max(0, D[i][j])));

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
		int strY = -1;

		// find start and end of alignment
		for (int i = 0; i != M.length; i++) {
			for (int j = 0; j != M[i].length; j++) {
				if (max <= M[i][j]) {
					max = M[i][j];
					x = i - 1;
					y = j - 1;
				}
			}
		}

		int score = max;
		String row0 = "";
		String row1 = "";
		int actScore = 0;
		char actx;
		char acty;

		// now cover the regions up until x and y - the gaps at the start and
		// the end
		for (int i = M[M.length - 1].length - 1; i > y + 1; i--) {
			row0 += "-";
			row1 += sequence2.getComp(i - 1);
		}
		for (int i = M.length - 1; i > x + 1; i--) {
			row0 += sequence1.getComp(i - 1);
			row1 += "-";
		}

		while (x >= 0 && y >= 0 && M[x + 1][y + 1] != 0) {

			actScore = M[x + 1][y + 1];
			actx = (Character) sequence1.getComp(x);
			acty = (Character) sequence2.getComp(y);
			strY = secStruct[y];

			if (actScore == M[x][y] + match(actx, acty, strY)) {
				row0 += actx;
				row1 += acty;
				y--;
				x--;
			} else if (actScore == D[x + 1][y + 1]) {
				while (D[x + 1][y + 1] == D[x + 1][y] + gapExtend[strY]
						&& y > 0) {
					row0 += "-";
					row1 += acty;
					y--;
					acty = (Character) sequence2.getComp(y);
				}
				row0 += "-";
				row1 += acty;
				y--;
			} else if (actScore == I[x + 1][y + 1]) {
				while (I[x + 1][y + 1] == I[x][y + 1] + gapExtend[strY]
						&& x > 0) {
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
				flip(row1.toCharArray()), 1.0d * score / (Gotoh.FACTOR * Gotoh.FACTOR));
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
