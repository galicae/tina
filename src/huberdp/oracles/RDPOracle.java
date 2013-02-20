/**
 * 
 */
package huberdp.oracles;

import huberdp.Oracle;
import huberdp.PartialAlignment;
import huberdp.RDPProblem;
import huberdp.scoring.RDPScoring;

import java.util.LinkedList;

import bioinfo.Sequence;
import bioinfo.alignment.Alignment;
import bioinfo.alignment.SequenceAlignment;
import bioinfo.alignment.Threading;
import bioinfo.alignment.gotoh.Gotoh;
import bioinfo.alignment.gotoh.LocalSequenceGotoh;
import bioinfo.proteins.PDBEntry;

/**
 * @author huberste
 * 
 */
public class RDPOracle extends LocalSequenceGotoh implements Oracle {

	private RDPScoring scoring;
	private PDBEntry template;
	private Sequence target;
	
	/**
	 * for local gotoh
	 */
	private boolean[][][] comefrom;

	public RDPOracle(double gapOpen, double gapExtend,
			double[][] scoringmatrix, RDPScoring scoring) {
		super(gapOpen, gapExtend, scoringmatrix);
		this.scoring = scoring;
	}

	public RDPOracle() {
		this(0.0, 0.0, null, null);
	}

	@Override
	protected void calculateMatrices() {
		Threading threading = new Threading(template, target, null, 0.0);
		for (int i = 1; i <= sequence1.length(); i++) {
			for (int j = 1; j <= sequence2.length(); j++) {
				D[i][j] = Math.max(M[i][j - 1] + gapOpen + gapExtend,
						D[i][j - 1] + gapExtend);
				I[i][j] = Math.max(M[i - 1][j] + gapOpen + gapExtend,
						I[i - 1][j] + gapExtend);
				M[i][j] = Math.max(
						M[i - 1][j - 1]
								+ (int) scoring.getScore(threading, i - 1,
										j - 1) * Gotoh.FACTOR,
						Math.max(I[i][j], Math.max(D[i][j], 0)));

			}
		}
	}
	
	/**
	 * @return Alignment of the two given Alignables
	 */
	@Override
	protected Alignment traceback() {
		
		int max = Integer.MIN_VALUE;
		int x = 0;
		int y = 0;

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

			if (actScore == M[x][y] + score(actx, acty)) { // come from topleft
				row0 += actx;
				row1 += acty;
				y--;
				x--;
			} else if (actScore == D[x + 1][y + 1]) { // come from left (deletion)
				while (D[x + 1][y + 1] == D[x + 1][y] + gapExtend && y > 0) {
					row0 += "-";
					row1 += acty;
					y--;
					acty = (Character)sequence2.getComp(y);
				}
				row0 += "-";
				row1 += acty;
				y--;
			} else if (actScore == I[x + 1][y + 1]) { // come from top (insertion)
				while (I[x + 1][y + 1] == I[x][y + 1] + gapExtend && x > 0) {
					row0 += actx;
					row1 += "-";
					x--;
					actx = (Character)sequence1.getComp(x);
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

		return new SequenceAlignment((Sequence)sequence1, (Sequence)sequence2,
				flip(row0.toCharArray()), flip(row1.toCharArray()), 1.0d
						* score / Gotoh.FACTOR);
	}

	@Override
	public LinkedList<PartialAlignment> findSimiliarSegments(
			RDPProblem problem, int m) {

		// allocate result
		LinkedList<PartialAlignment> results = new LinkedList<PartialAlignment>();

		// set Sequences for Oracle
		String[] rows = problem.getThreading().getRowsAsString();
		String template = "";
		String target = "";
		for (int i = problem.getProblemStart(); i < problem.getProblemEnd(); i++) {
			if (rows[0].charAt(i) != '-') {
				template += rows[0].charAt(i);
			}
			if (rows[1].charAt(i) != '-') {
				target += rows[1].charAt(i);
			}
		}

		rows = null; // GC

		// TODO make Sequence Alignment via dynamic programming with scoring
		// function
		SequenceAlignment alignment = null;

		results.add(new PartialAlignment(problem, alignment));

		return results;
	}

}
