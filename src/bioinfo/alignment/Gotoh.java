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
	
	int gap = 0;
	int[][] M;
	int[][] I;
	int[][] D;
	
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
	public Gotoh(int gap){
		this.gap = gap;
	}
	
	@Override
	public Alignment align(Alignable sequence1, Alignable sequence2) {
		M = new int[sequence1.length()][sequence2.length()];
		I = new int[sequence1.length()][sequence2.length()];
		D = new int[sequence1.length()][sequence2.length()];
		prepareMatrices();
		calculateMatrices();
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
		//TODO
		return null;
	}
	
	/**
	 * prepares matrices for global alignment
	 * 
	 */
	private void prepareMatrices(){
		//TODO
	}
	
	/**
	 * calculates matrices using scoring-function and gap-penalty
	 * 
	 */
	private void calculateMatrices(){
		//TODO
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
