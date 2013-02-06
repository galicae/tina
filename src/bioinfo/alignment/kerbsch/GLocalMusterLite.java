package bioinfo.alignment.kerbsch;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

import db.mysql.DBConnector;
import db.mysql.LocalConnection;
import bioinfo.Sequence;
import bioinfo.alignment.Alignable;
import bioinfo.alignment.Alignment;
import bioinfo.alignment.SequenceAlignment;
import bioinfo.alignment.gotoh.Gotoh;
import bioinfo.alignment.kerbsch.temp.SeqLibrary;
import bioinfo.proteins.DSSPEntry;
import bioinfo.proteins.DSSPFileReader;
import bioinfo.proteins.SecStructEight;

public class GLocalMusterLite extends Gotoh {

	private static final int INIT_VAL = Integer.MIN_VALUE / 2;

	// size for matrices
	private int xsize;
	private int ysize;

	// substitution matrices
	private int[][] hbMatrix = new int[26][26];
	private int[][] polMatrix = new int[26][26];
	private int[][] secStructMatrix = new int[26][26];
	private int[][] substMatrix = new int[26][26];

	// feature weights
	private final int hbWeight;
	private final int polWeight;
	private final int secStructWeight;
	private final int substWeight;

	private int[][] tempScore;
	private int score;
	private int[][] alignedRes;
	
	private char[] seq1;
	private char[] seq2;
	private char[] sec1;
	private char[] sec2;

	private HashMap<String, char[]> secStructSeqlib;

	// private LocalConnection con = new LocalConnection();
	// private DBConnector dbconnector = new DBConnector(con);

	public GLocalMusterLite(double gapOpen, double gapExtend,
			HashMap<String, char[]> sslib, int[][] hbMatrix, int[][] polMatrix,
			int[][] secStructMatrix, int[][] substMatrix, int hbWeight,
			int polWeight, int ssWeight, int substWeight) {
		super(gapOpen, gapExtend);

		this.gapOpen *= Gotoh.FACTOR;
		this.gapExtend *= Gotoh.FACTOR;

		this.hbMatrix = hbMatrix;
		this.polMatrix = polMatrix;
		this.secStructMatrix = secStructMatrix;
		this.substMatrix = substMatrix;
		this.hbWeight = hbWeight;
		this.polWeight = polWeight;
		this.secStructWeight = ssWeight;
		this.substWeight = substWeight;
		this.secStructSeqlib = sslib;
		
		seq1 = ((Sequence) sequence1).getSequence();
		seq2 = ((Sequence) sequence2).getSequence();
		sec1 = secStructSeqlib.get(sequence1.getID());
		sec2 = secStructSeqlib.get(sequence2.getID());
	}

	public GLocalMusterLite(double gapOpen, double gapExtend,
			HashMap<String, char[]> sslib, double[][] hbMatrix,
			double[][] polMatrix, double[][] secStructMatrix,
			double[][] substMatrix, double hbWeight, double polWeight,
			double ssWeight, double substWeight) {
		super(gapOpen, gapExtend);

		this.gapOpen *= Gotoh.FACTOR;
		this.gapExtend *= Gotoh.FACTOR;

		for (int i = 0; i < substMatrix.length; i++) {
			for (int j = 0; j < substMatrix[i].length; j++) {
				this.substMatrix[i][j] = (int) (Gotoh.FACTOR * substMatrix[i][j]);
				this.hbMatrix[i][j] = (int) (Gotoh.FACTOR * hbMatrix[i][j]);
				this.polMatrix[i][j] = (int) (Gotoh.FACTOR * polMatrix[i][j]);
				this.secStructMatrix[i][j] = (int) (Gotoh.FACTOR * secStructMatrix[i][j]);
			}
		}

		this.hbWeight = (int) (hbWeight * Gotoh.FACTOR);
		this.polWeight = (int) (polWeight * Gotoh.FACTOR);
		this.secStructWeight = (int) (ssWeight * Gotoh.FACTOR);
		this.substWeight = (int) (substWeight * Gotoh.FACTOR);

		// this.hbWeight = 1;
		// this.polWeight = 1;
		// this.secStructWeight = 2;
		// this.substWeight = 2;
		this.secStructSeqlib = sslib;
	}

	public SequenceAlignment align(Alignable sequence1, Alignable sequence2) {
		this.xsize = sequence1.length() + 1;
		this.ysize = sequence2.length() + 1;

		this.alignedRes = new int[2][];
		this.alignedRes[0] = new int[xsize - 1];
		this.alignedRes[1] = new int[ysize - 1];

		this.M = new int[xsize][ysize];
		this.I = new int[xsize][ysize];
		this.D = new int[xsize][ysize];
		tempScore = new int[xsize][ysize];

		
		//readin sequences
		this.sequence1 = (Sequence) sequence1;
		this.sequence2 = (Sequence) sequence2;
		
		seq1 = ((Sequence) sequence1).getSequence();
		seq2 = ((Sequence) sequence2).getSequence();
		sec1 = secStructSeqlib.get(sequence1.getID());
		sec2 = secStructSeqlib.get(sequence2.getID());
		
		int maxtemp = Integer.MAX_VALUE;
		//make tempscore
		for (int i = 1; i < xsize; i++) {
			for (int j = 1; j < ysize; j++) {
				tempScore[i][j] = score(hbMatrix, seq1[i - 1], seq2[j - 1])
						* hbWeight + score(polMatrix, seq1[i - 1], seq2[j - 1])
						* polWeight
						+ score(secStructMatrix, sec1[i - 1], sec2[j - 1])
						* secStructWeight
						+ score(substMatrix, seq1[i - 1], seq2[j - 1])
						* substWeight;
				if(tempScore[i][j] < maxtemp){
					maxtemp = tempScore[i][j];
				}
			}
		}
//		System.out.println((double)maxtemp/(Gotoh.FACTOR*Gotoh.FACTOR));
				
		doAlignment(1, M.length-1, 1, M[0].length-1);	
		return makeAlignment();
	}
	
	private void doAlignment(int xStart, int xEnd, int yStart, int yEnd){
		prepareMatrices(xStart, xEnd, yStart, yEnd);
		calculateMatrices(xStart, xEnd, yStart, yEnd);
		traceback(xStart, xEnd, yStart, yEnd);
	}

	private void prepareMatrices(int xStart, int xEnd, int yStart, int yEnd) {
		for (int i = xStart; i <= xEnd; i++) {
			D[i][yStart-1] = INIT_VAL;
			M[i][yStart-1] = 0;
			I[i][yStart-1] = 0;
		}

		for (int i = yStart; i <= yEnd; i++) {
			I[xStart-1][i] = INIT_VAL;
			M[xStart-1][i] = 0;
			D[xStart-1][i] = 0;
		}
		M[xStart-1][yStart-1] = 0;
		D[xStart-1][yStart-1] = INIT_VAL;
		I[xStart-1][yStart-1] = INIT_VAL;
	}

	private int score(int[][] matrix, char x, char y) {
		return matrix[x - 65][y - 65];
	}

	private void calculateMatrices(int xStart, int xEnd, int yStart, int yEnd) {
		// DSSPEntry template = dbconnector.getDSSP(sequence1.getID());
		// SecStructEight[] tSecStruct = template.getSecondaryStructure();
		//
		// DSSPEntry query = dbconnector.getDSSP(sequence2.getID());
		// SecStructEight[] qSecStruct = query.getSecondaryStructure();

		for (int i = xStart; i <= xEnd; i++) {
			for (int j = yStart; j <= yEnd; j++) {
				D[i][j] = Math.max(M[i][j - 1] + gapOpen + gapExtend,
						D[i][j - 1] + gapExtend);
				I[i][j] = Math.max(M[i - 1][j] + gapOpen + gapExtend,
						I[i - 1][j] + gapExtend);
				M[i][j] = Math.max(M[i - 1][j - 1] + tempScore[i][j],
						Math.max(I[i][j], Math.max(D[i][j], 0)));
			}
		}		
	}

	private SequenceAlignment makeAlignment() {
		int pointerX = 0;
		int pointerY = 0;
		List<Character> row1 = new ArrayList<Character>();
		List<Character> row2 = new ArrayList<Character>();

		while(pointerX < xsize - 1 || pointerY < ysize - 1){
			while(pointerX < xsize - 1){
				if (alignedRes[0][pointerX] == -1) {
					row2.add('-');
					row1.add((Character)sequence1.getComp(pointerX));
					pointerX++;
				} else {
					if (alignedRes[0][pointerX] == pointerY) {
						row1.add((Character)sequence1.getComp(pointerX));
						row2.add((Character)sequence2.getComp(pointerY));
						pointerX++;
						pointerY++;
					} else {
						break;
					}
				}
			}
			while(pointerY < ysize - 1){
				if (alignedRes[1][pointerY] == -1) {
					row1.add('-');
					row2.add((Character)sequence2.getComp(pointerY));
					pointerY++;
				} else {
					if (alignedRes[1][pointerY] == pointerX) {
						row1.add((Character)sequence1.getComp(pointerX));
						row2.add((Character)sequence2.getComp(pointerY));
						pointerX++;
						pointerY++;
					} else {
						break;
					}
				}
			}
		}
		
		//creating a chart array of aligned sequences
		char[] row_1 = new char[row1.size()];
		char[] row_2 = new char[row2.size()];
		for(int i = 0; i < row1.size(); i++){
			row_1[i] = row1.get(i);
		}
		for(int i = 0; i < row2.size(); i++){
			row_2[i] = row2.get(i);
		}
		
		return new SequenceAlignment((Sequence) sequence1,
				(Sequence) sequence2, row_1, row_2, 1.0d * score
				/ (Gotoh.FACTOR*Gotoh.FACTOR));
	}

	private int[] findMaxScore(int fromX, int toX, int fromY, int toY) {
		int max = INIT_VAL;
		int maxX = 0;
		int maxY = 0;
		for (int i = fromX; i <= toX; i++) {
			for (int j = fromY; j <= toY; j++) {
				if (M[i][j] > max) {
					max = M[i][j];
					maxX = i;
					maxY = j;
				}
			}
		}
		return new int[] { maxX, maxY };
	}

	private void traceback(int xStart, int xEnd, int yStart, int yEnd) {

		int[] maxScore = findMaxScore(xStart, xEnd, yStart, yEnd);
		score += M[maxScore[0]][maxScore[1]];
		int x = maxScore[0];
		int y = maxScore[1];

		int actScore;
		
		while (x >= xStart && y >= yStart && M[x][y] != 0) {

			actScore = M[x][y];

			if (actScore == M[x - 1][y - 1] + tempScore[x][y]) {
				alignedRes[0][x - 1] = y-1;
				alignedRes[1][y - 1] = x-1;
				y--;
				x--;
			} else if (actScore == D[x][y]) {
				while (D[x][y] == D[x][y - 1] + gapExtend && y > 0) {
					alignedRes[1][y - 1] = -1;
					y--;
				}
				alignedRes[1][y - 1] = -1;
				y--;
			} else if (actScore == I[x][y]) {
				while (I[x][y] == I[x - 1][y] + gapExtend && x > 0) {
					alignedRes[0][x - 1] = -1;
					x--;
				}
				alignedRes[0][x - 1] = -1;
				x--;
			}
		}

		// search for local alignments in area after actual local alignment
		if (maxScore[0]+1 < xEnd && maxScore[1]+1 < yEnd) {
			doAlignment(maxScore[0]+1, xEnd, maxScore[1]+1, yEnd);
		} else {
			for (int i = maxScore[0]+1; i <= xEnd; i++) {
				alignedRes[0][i - 1] = -1;
			}
			for (int i = maxScore[1]+1; i <= yEnd; i++) {
				alignedRes[1][i - 1] = -1;
			}
		}

		// search for local alignments in area before actual local alignment
		if (xStart < x && yStart < y) {
			doAlignment(xStart, x, yStart, y);
		} else {
			for (int i = x; i >= xStart; i--) {
				alignedRes[0][i - 1] = -1;
			}
			for (int i = y; i >= yStart; i--) {
				alignedRes[1][i - 1] = -1;
			}
		}
	}
}
