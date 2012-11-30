package bioinfo.alignment;

import java.util.ArrayList;
import java.util.List;

import bioinfo.Sequence;
import bioinfo.proteins.AminoAcid;
import bioinfo.proteins.AtomType;
import bioinfo.proteins.PDBEntry;

/**
 * 
 * @author andreseitz
 * Freeshift Alignment of two Sqeuences
 */
public class EuklidPDBEntryGotoh extends Gotoh{

	private static final int INIT_VAL = Integer.MIN_VALUE/2;

	/**
	 * 
	 * @param gapOpen
	 * @param gapExtend
	 */
	public EuklidPDBEntryGotoh(int gapOpen, int gapExtend) {
		super(gapOpen*1000, gapExtend*1000);
	}
	
	@Override
	public StructureAlignment align(Alignable sequence1, Alignable sequence2) {
		this.M = new int[sequence1.length()+1][sequence2.length()+1];
		this.I = new int[sequence1.length()+1][sequence2.length()+1];
		this.D = new int[sequence1.length()+1][sequence2.length()+1];
		this.sequence1 = (PDBEntry)sequence1;
		this.sequence2 = (PDBEntry)sequence2;
		prepareMatrices();
		calculateMatrices();
		return (StructureAlignment)traceback();
	}

	@Override
	public boolean check(Alignment alignment) {
		StructureAlignment ali = (StructureAlignment)alignment;
		int score = 0;
		AminoAcid[] row0 = ali.getRow(0);
		AminoAcid[] row1 = ali.getRow(1);
		for(int i = 0; i < row0.length; i++){
			if(row0[i] == null || row1[i] == null){
				score += gapOpen;
				while(i < row0.length && (row0[i] == null || row1[i] == null)){
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
				M[i][j] = Math.max(M[i-1][j-1]+score((AminoAcid)sequence1.getComp(i-1),(AminoAcid)sequence2.getComp(j-1)),
							Math.max(I[i][j],
							D[i][j]));
				
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
			if(M[i][M[i].length-1] >= max){
				max = M[i][M[i].length-1];
				x = i-1;
				y = M[i].length-2;
			}
		}
		for(int i = (M[M.length-1].length-1); i >= 0; i--){
			if(M[(M.length-1)][i] > max){
				max = M[(M.length-1)][i];
				y = i-1;
				x = M.length-2;
			}
		}
		
		int score = max;
		List<AminoAcid> row0 = new ArrayList<AminoAcid>();
		List<AminoAcid> row1 = new ArrayList<AminoAcid>();
		int actScore = 0;
		AminoAcid actx;
		AminoAcid acty;
		
		for(int i = M[M.length-1].length-1; i > y+1; i--){
			row0.add(null);
			row1.add((AminoAcid)sequence2.getComp(i-1));
		}
		for(int i = M.length-1; i > x+1; i--){
			row0.add((AminoAcid)sequence1.getComp(i-1));
			row1.add(null);
		}
		
		while(x >= 0 && y >= 0){
			
			actScore = M[x+1][y+1];
			actx = (AminoAcid)sequence1.getComp(x);
			acty = (AminoAcid)sequence2.getComp(y);
			
			if(actScore == M[x][y]+score(actx, acty)){
				row0.add(actx);
				row1.add(acty);
				y--;
				x--;
			} else if (actScore == D[x+1][y+1]){
				while(D[x+1][y+1] == D[x+1][y]+gapExtend && y > 0){
					row0.add(null);
					row1.add(acty);
					y--;
					acty = (AminoAcid)sequence2.getComp(y);
				}
				row0.add(null);
				row1.add(acty);
				y--;
			} else if (actScore == I[x+1][y+1]){
				while(I[x+1][y+1] == I[x][y+1]+gapExtend && x > 0){
					row0.add(actx);
					row1.add(null);
					x--;
					actx = (AminoAcid)sequence1.getComp(x);
				}
				row0.add(actx);
				row1.add(null);
				x--;
			}	
		}
		for(int i = y+1; i > 0; i--){
			row0.add(null);
			row1.add((AminoAcid)sequence2.getComp(i-1));
		}
		for(int i = x+1; i > 0; i--){
			row0.add((AminoAcid)sequence1.getComp(i-1));
			row1.add(null);
		}
		
		return new StructureAlignment((PDBEntry)sequence1, (PDBEntry)sequence2, flip(row0.toArray(new AminoAcid[row0.size()])), flip(row1.toArray(new AminoAcid[row1.size()])), score);
	}

	/**
	 * @param two components of Alignable implementing equals
	 * @return score between two components of Alignable
	 */
	private int score(AminoAcid x, AminoAcid y) {
		double[] p = x.getAtomByType(AtomType.CA).getPosition();
		double[] q = y.getAtomByType(AtomType.CA).getPosition();
		return (int)(Math.sqrt(((p[0]-q[0])*(p[0]-q[0]))+((p[1]-q[1])*(p[1]-q[1]))+((p[2]-q[2])*(p[2]-q[2])))*1000);
	}
	
	private AminoAcid[] flip(AminoAcid[] in){
		AminoAcid[] out = new AminoAcid[in.length];
		for(int i = in.length-1; i >= 0; i --){
			out[out.length-1-i] = in[i];
		}
		return out;
	}

}

