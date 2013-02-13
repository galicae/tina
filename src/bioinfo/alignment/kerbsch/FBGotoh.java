package bioinfo.alignment.kerbsch;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;

import bioinfo.Sequence;
import bioinfo.alignment.Alignable;
import bioinfo.alignment.kerbsch.temp.LocalCore;
import bioinfo.alignment.kerbsch.temp.LocalMatch;

public class FBGotoh {
	private static final int INIT_VAL = Integer.MIN_VALUE / 2;
	private final int FACTOR = 100;
	private int gapOpen;
	private int gapExtend;
	private int[][] M;
	private int[][] I;
	private int[][] D;
	
	private int[][] substMatrix;
	private int[][] revM;
	private int[][] revD;
	private int[][] revI;
	private int[][] hybridM;
	
	// size for matrices
	private int xsize;
	private int ysize;

	// for traceback
	private int[][] tempScore;
	private int[][] alignedRes;

	private char[] seq1;
	private char[] seq2;

	
	//other stuff
	private BufferedWriter out;
	
	
	public FBGotoh(double gapOpen, double gapExtend, int[][] scoringmatrix, BufferedWriter out) {
		this.gapOpen = (int) (gapOpen * FACTOR);
		this.gapExtend = (int) (gapExtend * FACTOR);
		this.substMatrix = scoringmatrix;
		this.out = out;
	}

	public FBGotoh(double gapOpen, double gapExtend, double[][] scoringmatrix, BufferedWriter out) {
		this.gapOpen = (int) (gapOpen * FACTOR);
		this.gapExtend = (int) (gapExtend * FACTOR);
		
		this.substMatrix = new int[scoringmatrix.length][scoringmatrix[0].length];
		for (int i = 0; i != scoringmatrix.length; i++) {
			for (int j = 0; j != scoringmatrix[0].length; j++) {
				this.substMatrix[i][j] = (int) (FACTOR * scoringmatrix[i][j]);
			}
		}
		
		this.out = out;
	}

	public int[][] align(Alignable sequence1, Alignable sequence2) {
		this.xsize = sequence1.length() + 1;
		this.ysize = sequence2.length() + 1;
		
		this.alignedRes = new int[2][];
		this.alignedRes[0] = new int[xsize - 1];
		this.alignedRes[1] = new int[ysize - 1];

		this.M = new int[xsize][ysize];
		this.I = new int[xsize][ysize];
		this.D = new int[xsize][ysize];
		this.revM = new int[xsize][ysize];
		this.revI = new int[xsize][ysize];
		this.revD = new int[xsize][ysize];
		this.hybridM = new int[xsize][ysize];
		tempScore = new int[xsize][ysize];
		
		// readin sequences
		seq1 = ((Sequence)sequence1).getSequence();
		seq2 = ((Sequence)sequence2).getSequence();

		// make tempscore
		for (int i = 1; i < xsize; i++) {
			for (int j = 1; j < ysize; j++) {
				tempScore[i][j] = score(seq1[i - 1], seq2[j - 1]);
			}
		}
		prepareMatrices();
		calculateMatrices();
		findLocals();
		return this.hybridM;
	}

	private void prepareMatrices() {
		for (int i = 0; i < xsize; i++) {
			D[i][0] = INIT_VAL;
			revD[i][0] = INIT_VAL;
		}

		for (int i = 0; i < ysize; i++) {
			I[0][i] = INIT_VAL;
			revI[0][i] = INIT_VAL;
		}
	}

	private int score(char x, char y) {
		return substMatrix[x - 65][y - 65];
	}

	private void calculateMatrices() {

		for (int i = 1; i < xsize; i++) {
			for (int j = 1; j < ysize; j++) {
				D[i][j] = Math.max(M[i][j - 1] + gapOpen + gapExtend,
						D[i][j - 1] + gapExtend);
				revD[i][j] = Math.max(revM[i][j - 1] + gapOpen + gapExtend,
						revD[i][j - 1] + gapExtend);

				I[i][j] = Math.max(M[i - 1][j] + gapOpen + gapExtend,
						I[i - 1][j] + gapExtend);
				revI[i][j] = Math.max(revM[i - 1][j] + gapOpen + gapExtend,
						revI[i - 1][j] + gapExtend);

				revM[i][j] = Math.max(revM[i - 1][j - 1]
						+ tempScore[xsize - i][ysize - j],
						Math.max(revI[i][j], Math.max(revD[i][j], 0)));
				M[i][j] = Math.max(M[i - 1][j - 1] + tempScore[i][j],
						Math.max(I[i][j], Math.max(D[i][j], 0)));
			}
		}
		
		for (int i = 1; i < xsize; i++) {
			for (int j = 1; j < ysize; j++) {
				hybridM[i][j] = M[i][j] + revM[xsize - i][ysize - j]; 
			}
		}		
	}
	
	private int calcScoreAndLength(int startIndex, int endIndex, List<LocalMatch> localmatches){
		int[] result = new int[2];
		int score = 0;
		int length = endIndex-startIndex+1;
		
//		for (int i = startIndex; i <= endIndex; i++) {
//			
//		}
		
		result[0] = score;
		result[1] = length;
		return length;
	}
	
	private List<LocalCore> traceback(List<LocalMatch> localmatches){
		int x,y,revX,revY;
		int actScore;
		List<LocalCore> result = new ArrayList<LocalCore>();
		ArrayList<int[]> lastcoreCoords;
		
		int[][] used = new int[xsize][ysize];

		for (LocalMatch lm : localmatches) {
			
			x = lm.getCoords()[0];
			y = lm.getCoords()[1];	
			revX = xsize - x;
			revY = ysize - y;
			
			//do reverse traceback
			if(used[x][y] == -1){
				continue;
			}else{
				result.add(new LocalCore(lm.getScore()));
				lastcoreCoords = result.get(result.size()-1).getCoords();
				lastcoreCoords.add(new int[]{xsize-revX,ysize-revY});
				used[xsize-revX][ysize-revY] = -1;
				
				while(revX > 0 && revY > 0 && revM[revX][revY] != 0){
					actScore = revM[revX][revY];
					
					if (actScore == revM[revX - 1][revY - 1] + tempScore[xsize-revX][ysize-revY]) {
						revY--;
						revX--;
						if(used[xsize-revX][ysize-revY] != -1){
							lastcoreCoords.add(new int[]{xsize-revX,ysize-revY});
							used[xsize-revX][ysize-revY] = -1;
						}else{
							break;
						}
						
					} else if (actScore == revD[revX][revY]) {
						while (revY > 0 && revD[revX][revY] == revD[revX][revY - 1] + gapExtend) {
							revY--;
							if(used[xsize-revX][ysize-revY] != -1){
								lastcoreCoords.add(new int[]{xsize-revX,ysize-revY});
								used[xsize-revX][ysize-revY] = -1;
							} else {
								break;
							}
						}
						revY--;
						if(used[xsize-revX][ysize-revY] != -1){
							lastcoreCoords.add(new int[]{xsize-revX,ysize-revY});
							used[xsize-revX][ysize-revY] = -1;	
						} else {
							break;
						}
						
					} else if (actScore == revI[revX][revY]) {
						while (revX > 0 && revI[revX][revY] == revI[revX - 1][revY] + gapExtend) {					
							revX--;
							if(used[xsize-revX][ysize-revY] != -1){
								lastcoreCoords.add(new int[]{xsize-revX,ysize-revY});
								used[xsize-revX][ysize-revY] = -1;	
							} else {
								break;
							}
						}
						revX--;					
						if(used[xsize-revX][ysize-revY] != -1){
							lastcoreCoords.add(new int[]{xsize-revX,ysize-revY});
							used[xsize-revX][ysize-revY] = -1;
						} else {
							break;
						}
					}
				}
				
				//do normal traceback
				while(x > 0 && y > 0 && M[x][y] != 0){
					actScore = M[x][y];
					
					if (actScore == M[x - 1][y - 1] + tempScore[x][y]) {
						y--;
						x--;
						if(used[x][y] != -1){
							lastcoreCoords.add(new int[]{x,y});
							used[x][y] = -1;
						}else{
							break;
						}
						
					} else if (actScore == D[x][y]) {
						while (D[x][y] == D[x][y - 1] + gapExtend && y > 0) {
							y--;
							if(used[x][y] != -1){
								lastcoreCoords.add(new int[]{x,y});
								used[x][y] = -1;
							} else {
								break;
							}
						}
						y--;
						if(used[x][y] != -1){
							lastcoreCoords.add(new int[]{x,y});
							used[x][y] = -1;	
						} else {
							break;
						}
						
					} else if (actScore == I[x][y]) {
						while (I[x][y] == I[x - 1][y] + gapExtend && x > 0) {					
							x--;
							if(used[x][y] != -1){
								lastcoreCoords.add(new int[]{x,y});
								used[x][y] = -1;	
							} else {
								break;
							}
						}
						x--;					
						if(used[x][y] != -1){
							lastcoreCoords.add(new int[]{x,y});
							used[x][y] = -1;
						} else {
							break;
						}
					}
				}
			}
		}
		return result;
	}
	
	private class sortScore implements Comparator<LocalMatch>{

		@Override
		public int compare(LocalMatch arg0, LocalMatch arg1) {
			return arg1.getScore() - arg0.getScore();
		}
		
	}
	
	private void findLocals(){
		int lengthCutOff = 3;
		int scoreCutOff = (int) (10.0 * FACTOR);
		List<LocalMatch> localmatches = new ArrayList<LocalMatch>();
		
		
		
		//read all scores with their coordinates
		for (int i = xsize-1; i > 0; i--) {
			for (int j = ysize-1; j > 0; j--) {
				if(hybridM[i][j] >= scoreCutOff){
					localmatches.add(new LocalMatch(hybridM[i][j],new int[]{i,j}));
				}
			}
		}

		//sort scorelist
		Collections.sort(localmatches, new sortScore());

		List<LocalCore> localcores = traceback(localmatches);

		//sort scorelist
//		Collections.sort(localmatches, new sortScore());
				
		try {
//			int lastscore = localmatches.get(0).getScore();
//			int actscore;
//			int startIndex = 0;
//			int endIndex = 0;
//			int result;
			for (LocalCore ue : localcores) {
				if(ue.getCoords().size() >= lengthCutOff){
					out.write(ue.getScore() + ": "+ ue.getCoords().size()+"\n");
				}
//				actscore = ue.getScore();
//				if(actscore == lastscore){
//					endIndex++;
//				} else {	
//					result = calcScoreAndLength(startIndex,endIndex,localmatches);
//					if(result >= lengthCutOff){
//						out.write(lastscore + ": "+ result+"\n");
//					}
//					lastscore = actscore;
//					endIndex++;
//					startIndex = endIndex;
//				}
			}
		} catch (IOException e) {
			System.out.println("failure");
		}
		
	}
	
	public void printHM(){
		for (int i = 0; i < xsize; i++) {
			for (int j = 0; j < ysize; j++) {
				System.out.print(hybridM[i][j]+"\t");
			}
			System.out.println();
		}
	}
	
	public void printM(){
		for (int i = 0; i < xsize; i++) {
			for (int j = 0; j < ysize; j++) {
				System.out.print(M[i][j]+"\t");
			}
			System.out.println();
		}
	}
}
