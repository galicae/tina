package bioinfo.alignment;

import bioinfo.Sequence;

/**
 * 
 * @author andreseitz
 * Local Alignment of two Sequences
 */
public class LocalSequenceGotoh extends Gotoh{

	private static final int INIT_VAL = Integer.MIN_VALUE/2;
	private int[][] scoringmatrix;
	Sequence sequence1;
	Sequence sequence2;
	
	/**
	 * 
	 * @param gapOpen
	 * @param gapExtend
	 * @param scoringmatrix 26x26 matrix containing all scoring values plus some empty lines due to faster access
	 */
	public LocalSequenceGotoh(int gapOpen, int gapExtend, int[][] scoringmatrix) {
		super(gapOpen, gapExtend);
		this.scoringmatrix = scoringmatrix;
	}
	
	@Override
	public SequenceAlignment align(Alignable sequence1, Alignable sequence2) {
		this.M = new int[sequence1.length()+1][sequence2.length()+1];
		this.I = new int[sequence1.length()+1][sequence2.length()+1];
		this.D = new int[sequence1.length()+1][sequence2.length()+1];
		this.sequence1 = (Sequence)sequence1;
		this.sequence2 = (Sequence)sequence2;
		prepareMatrices();
		calculateMatrices();
		return (SequenceAlignment)traceback();
	}

	@Override
	public boolean check(Alignment alignment) {
		SequenceAlignment ali = (SequenceAlignment)alignment;
		int score = 0;
		char[] row0 = ali.getRow(0);
		char[] row1 = ali.getRow(1);
		for(int i = 0; i < row0.length; i++){
			if(row0[i] == '-' || row1[i] == '-'){
				score += gapOpen;
				while(i < row0.length && (row0[i] == '-' || row1[i] == '-')){
					score += this.gapExtend;
					i++;
				}
				i--;
			}else {
				score += score(row0[i], row1[i]);	
			}
		}
		if(score == ali.getScore()){
			return true;
		}else{
			return false;
		}
	}
	
	/**
	 * prepares matrices for global alignment
	 * 
	 */
	private void prepareMatrices() {
		for(int i = 1; i <= sequence1.length(); i++){ //old: hSeq
			//I[i][0] = 0; //old: vGap
			D[i][0] = INIT_VAL; //old: hGap
			//M[i][0] = 0;
		}
		
		for(int i = 1; i <= sequence2.length(); i++){ //old: vSeq
			//D[0][i] = 0;
			//M[0][i] = 0;
			I[0][i] = INIT_VAL;
		}
		D[0][0] = INIT_VAL;
		I[0][0] = INIT_VAL;
		//M[0][0] = 0;		
	}
	
	/**
	 * calculates matrices using scoring-function and gap-penalty
	 * 
	 */
	private void calculateMatrices() {
		for(int i = 1; i <= sequence1.length(); i++){
			for(int j = 1; j <= sequence2.length(); j++){
				D[i][j] = Math.max(M[i][j-1]+gapOpen, D[i][j-1]+gapExtend);
				I[i][j] = Math.max(M[i-1][j]+gapOpen,I[i-1][j]+gapExtend);
				M[i][j] = Math.max(M[i-1][j-1]+score(sequence1.getComp(i-1),sequence2.getComp(j-1)),
							Math.max(I[i][j],
							Math.max(D[i][j],0)));
				
			}
		}		
	}
	
	/**
	 * Override this method in extensions!
	 * @return Alignment of the two given Alignables
	 */
	private Alignment traceback() {
		
		int max = INIT_VAL;
		int x =	0;
		int y = 0;
		
		for(int i = 0; i != M.length; i++){
			for(int j = 0; j != M[i].length; j++){
				if(max <= M[i][j]){
					max = M[i][j];
					x = i-1;
					y = j-1;
				}
			}
		}
		
		int score = max;
		String row0 = "";
		String row1 = "";
		int actScore = 0;
		char actx;
		char acty;
		
		for(int i = M[M.length-1].length-1; i > y+1; i--){
			row0 += "-";
			row1 += sequence2.getComp(i-1);
		}
		for(int i = M.length-1; i > x+1; i--){
			row0 += sequence1.getComp(i-1);
			row1 += "-";
		}
		
		while(x >= 0 && y >= 0 && M[x+1][y+1] != 0){
			
			actScore = M[x+1][y+1];
			actx = sequence1.getComp(x);
			acty = sequence2.getComp(y);
			
			if(actScore == M[x][y]+score(actx, acty)){
				row0 += actx;
				row1 += acty;
				y--;
				x--;
			} else if (actScore == D[x+1][y+1]){
				while(D[x+1][y+1] == D[x+1][y]+gapExtend && y > 0){
					row0 += "-";
					row1 += acty;
					y--;
					acty = sequence2.getComp(y);
				}
				row0 += "-";
				row1 += acty;
				y--;
			} else if (actScore == I[x+1][y+1]){
				while(I[x+1][y+1] == I[x][y+1]+gapExtend && x > 0){
					row0 += actx;
					row1 += "-";
					x--;
					actx = sequence1.getComp(x);
				}
				row0 += actx;
				row1 += "-";
				x--;
			}	
		}
		for(int i = y+1; i > 0; i--){
			row0 += "-";
			row1 += sequence2.getComp(i-1);
		}
		for(int i = x+1; i > 0; i--){
			row0 += sequence1.getComp(i-1);
			row1 += "-";
		}
		
		return new SequenceAlignment(sequence1, sequence2, flip(row0.toCharArray()), flip(row1.toCharArray()), score);
	}

	/**
	 * @param two components of Alignable implementing equals
	 * @return score between two components of Alignable
	 */
	private int score(char x, char y) {
		return scoringmatrix[x-65][y-65];
	}
	
	private char[] flip(char[] in){
		char[] out = new char[in.length];
		for(int i = in.length-1; i >= 0; i --){
			out[out.length-1-i] = in[i];
		}
		return out;
	}
}
