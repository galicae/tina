package bioinfo.alignment;

import java.io.BufferedWriter;
import java.io.Writer;

/**
 * 
 * @author gobi4
 * abstract class Gotoh defines interface for alignments
 * if no method is overridden in extending class it performs an 
 * gotoh alignment on Einheits-Matrix with global matrix-preparation
 * and traceback
 * 
 * for other alignmenttypes the following methods have to be overridden:
 * 		+
 *
 */

public abstract class Gotoh implements Aligner {
	
	private int gapOpen = 0;
	private int gapExtend = 0;
	private int[][] M;
	private int[][] I;
	private int[][] D;
	
	/**
	 * Constructor initialising scoring by 
	 * Einheitsscoring (gap:0, subst:0, match:1) 
	 * if scoring function is not overridden
	 */
	public Gotoh(){
		
	}
	
	/**
	 * Constructor initialising Gotoh with
	 * certain gap-penalty
	 * @param gap integer defining gap-costs
	 */
	public Gotoh(int gapOpen, int gapExtend){
		this.gapOpen = gapOpen;
		this.gapExtend = gapExtend;
	}
	
	@Override
	public Alignment align(Alignable sequence1, Alignable sequence2) {
		M = new int[sequence1.length()][sequence2.length()];
		I = new int[sequence1.length()][sequence2.length()];
		D = new int[sequence1.length()][sequence2.length()];
		prepareMatrices(sequence1,sequence2);
		calculateMatrices(sequence1,sequence2);
		return traceback();
	}

	@Override
	public boolean check(Alignment arg) {
		// TODO Auto-generated method stub
		return false;
	}
	
	
	/**
	 * 
	 * @return Alignment of the two given Alignables
	 */
	private Alignment traceback(){
		
		return null;
	}
	
	/**
	 * prepares matrices for global alignment
	 * 
	 */
	private void prepareMatrices(Alignable sequence1, Alignable sequence2){
		for(int i = 1; i <= sequence1.length(); i++){ //old: hSeq
			I[i][0] = 0; //old: vGap
			D[i][0] = Integer.MIN_VALUE; //old: hGap
			if(i==1){
				M[i][0] = M[i-1][0]+gapOpen;
			}else{
				M[i][0] = M[i-1][0]+gapExtend;
			}
		}
		
		for(int i = 1; i <= sequence2.length(); i++){ //old: vSeq
			D[0][i] = 0;
			if(i==1){
				M[0][i] = M[0][i-1]+gapOpen;
			}else{
				M[0][i] = M[0][i-1]+gapExtend;
			}
			I[0][i] = Integer.MIN_VALUE;
		}
		D[0][0] = Integer.MIN_VALUE;
		I[0][0] = Integer.MIN_VALUE;
		M[0][0] = 0;
	}
	
	/**
	 * calculates matrices using scoring-function and gap-penalty
	 * 
	 */
	private void calculateMatrices(Alignable sequence1, Alignable sequence2){
		for(int i = 1; i <= sequence1.length(); i++){
			for(int j = 1; j <= sequence2.length(); j++){
				D[i][j] = Math.max(M[i][j-1]+gapOpen, D[i][j-1]+gapExtend);
				I[i][j] = Math.max(M[i-1][j]+gapOpen,I[i-1][j]+gapExtend);
				M[i][j] = Math.max(M[i-1][j-1]+score(sequence2.getComp(i-1),sequence1.getComp(j-1)),
							Math.max(I[i][j],
							D[i][j]));
				
			}
		}	
	}
	
	/**
	 * @param two components of Alignable implementing equals
	 * @return score between two components of Alignable
	 */
	private int score(Object x, Object y){
		if(x.equals(y)){
			return 1;
		}else{
			return 0;
		}
	}
	
	/**
	 * @param out Writer for example BufferedWriter to file or to System.out
	 * streams all matrices as tab-separated content
	 */
	public void streamMatricesAsTxt(Writer out){
		//TODO
	}
	
	/**
	 * @param out Writer for example BufferedWriter to file or to System.out
	 * streams all matrices as html-formatted content
	 */
	public void streamMatricesAsHtml(Writer out){
		//TODO
	}
	
	
	
	
}
