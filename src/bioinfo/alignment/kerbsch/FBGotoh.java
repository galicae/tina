package bioinfo.alignment.kerbsch;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;

public class FBGotoh {
	private static final int INIT_VAL = Integer.MIN_VALUE / 2;
	private final double FACTOR = 100.0;
	private final double LAMBDA;
	private final double DBLENGTH;
	private int lengthCutOff = 3;
	private double scoreCutOff = 0.0;
	private int gapOpen;
	private int gapExtend;
	private int[][] substMatrix;

	private int[][] M;
	private int[][] I;
	private int[][] D;
	private int[][] revM;
	private int[][] revD;
	private int[][] revI;
	private int[][] hybridM;

	// size for matrices
	private int xsize;
	private int ysize;

	// for traceback
	private int[][] tempScore;

	private char[] seq1;
	private char[] seq2;

	BufferedWriter out = new BufferedWriter(new FileWriter(
			"../Evalues_fbaligns/dayhoff_local.evalues"));

	public FBGotoh(double gapOpen, double gapExtend, double lambda,
			int dblength, int[][] scoringmatrix) throws IOException {
		this.gapOpen = (int) (gapOpen * FACTOR);
		this.gapExtend = (int) (gapExtend * FACTOR);
		this.substMatrix = scoringmatrix;
		DBLENGTH = dblength;
		LAMBDA = lambda;
	}

	public FBGotoh(double gapOpen, double gapExtend, double lambda,
			int dblength, double[][] scoringmatrix) throws IOException {
		this.gapOpen = (int) (gapOpen * FACTOR);
		this.gapExtend = (int) (gapExtend * FACTOR);
		DBLENGTH = dblength;
		LAMBDA = lambda;

		this.substMatrix = new int[scoringmatrix.length][scoringmatrix[0].length];
		for (int i = 0; i != scoringmatrix.length; i++) {
			for (int j = 0; j != scoringmatrix[0].length; j++) {
				this.substMatrix[i][j] = (int) (FACTOR * scoringmatrix[i][j]);
			}
		}
	}

	public List<Locals> align(char[] seq1, char[] seq2) {
		this.xsize = seq1.length + 1;
		this.ysize = seq2.length + 1;

		this.seq1 = seq1;
		this.seq2 = seq2;

		this.M = new int[xsize][ysize];
		this.I = new int[xsize][ysize];
		this.D = new int[xsize][ysize];
		this.revM = new int[xsize][ysize];
		this.revI = new int[xsize][ysize];
		this.revD = new int[xsize][ysize];
		this.hybridM = new int[xsize][ysize];
		tempScore = new int[xsize][ysize];

		// make tempscore
		for (int i = 1; i < xsize; i++) {
			for (int j = 1; j < ysize; j++) {
				tempScore[i][j] = score(seq1[i - 1], seq2[j - 1]);
			}
		}
		prepareMatrices();
		return calculateMatrices();
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

	private List<Locals> calculateMatrices() {

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

				M[i][j] = Math.max(M[i - 1][j - 1] + tempScore[i][j],
						Math.max(I[i][j], Math.max(D[i][j], 0)));
				revM[i][j] = Math.max(revM[i - 1][j - 1]
						+ tempScore[xsize - i][ysize - j],
						Math.max(revI[i][j], Math.max(revD[i][j], 0)));
			}
		}

		// main step, calculate hybridM, find locals, do traceback
		List<LocalMatch> scorelist = new ArrayList<LocalMatch>();
		for (int i = 1; i < xsize; i++) {
			for (int j = 1; j < ysize; j++) {
				hybridM[i][j] = M[i - 1][j - 1] + tempScore[i][j]
						+ revM[xsize - i - 1][ysize - j - 1];
				if (hybridM[i][j] >= scoreCutOff) {
					scorelist.add(new LocalMatch(hybridM[i][j], new int[] { i,
							j }));

				}
				// System.out.print(hybridM[i][j]+"\t");
			}
			// System.out.println();
		}

		// sort scorelist
		Collections.sort(scorelist, new sortScore());
		System.out.print("Maxscore: " + (scorelist.get(0).getScore() / 100.0)
				+ "\t");

		List<Locals> locals = traceback(scorelist);

		return locals;
	}

	private Locals mergeAlignment(List<int[]> coordsBW, List<int[]> coordsFW) {
		List<int[]> merged = new ArrayList<int[]>();

		// add forward alignment
		int lastindex = coordsFW.size() - 1;
		// while (lastindex > 0
		// && (coordsFW.get(lastindex - 1)[1] == coordsFW.get(lastindex)[1] ||
		// coordsFW
		// .get(lastindex - 1)[0] == coordsFW.get(lastindex)[0])) {
		// lastindex--;
		// }

		for (int i = lastindex; i >= 0; i--) {
			merged.add(coordsFW.get(i));
		}

		// add backward alignment

		lastindex = coordsBW.size() - 1;
		// while (lastindex > 0
		// && (coordsBW.get(lastindex - 1)[1] == coordsBW.get(lastindex)[1] ||
		// coordsBW
		// .get(lastindex - 1)[0] == coordsBW.get(lastindex)[0])) {
		// lastindex--;
		// }

		for (int i = 0; i <= lastindex; i++) {
			merged.add(coordsBW.get(i));
		}

		double evalue = Double.NEGATIVE_INFINITY;
		double newscore = calcScore(merged) / FACTOR;
		if (merged.size() >= lengthCutOff && newscore >= scoreCutOff) {
			evalue = Math.log(1.0d * DBLENGTH * seq1.length * seq2.length
					* Math.exp(-LAMBDA * newscore));
			try {
				out.append(((Math.round(evalue * 100.0)) / 100.0) + "\t"
						+ newscore + "\t" + merged.size() + "\n");
				// out.append(newscore + "\n");
				// out.append("\n");
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}
		return new Locals(newscore, evalue, merged);
	}

	private int calcScore(List<int[]> match) {
		int score = 0;
		int start = 0;
		int end = match.size() - 1;
		while (start < end
				&& (match.get(start + 1)[0] == match.get(start)[0] || match
						.get(start + 1)[1] == match.get(start)[1])) {
			start++;
		}

		while (end > 1
				&& (match.get(end - 1)[0] == match.get(end)[0] || match
						.get(end - 1)[1] == match.get(end)[1])) {
			end--;
		}

		for (int i = start; i <= end; i++) {
			score += tempScore[match.get(i)[0]][match.get(i)[1]];
			if (i < end - 1
					&& (match.get(i + 1)[0] == match.get(i)[0] || match
							.get(i + 1)[1] == match.get(i)[1])) {
				score += gapOpen;
				while (i < end - 1
						&& (match.get(i + 1)[0] == match.get(i)[0] || match
								.get(i + 1)[1] == match.get(i)[1])) {
					i++;
					score += gapExtend;
				}
				i--;
			}
		}
		return score;
	}

	private List<Locals> traceback(List<LocalMatch> localmatches) {
		int x, y, revX, revY;
		int actScore;
		double maxscore = Double.NEGATIVE_INFINITY;

		List<Locals> result = new ArrayList<Locals>();
		ArrayList<int[]> coordsFW = new ArrayList<int[]>();
		ArrayList<int[]> coordsBW = new ArrayList<int[]>();

		int[][] used = new int[xsize][ysize];

		for (LocalMatch lm : localmatches) {
			coordsBW = new ArrayList<int[]>();
			coordsFW = new ArrayList<int[]>();

			x = lm.getCoords()[0];
			y = lm.getCoords()[1];

			// do reverse traceback
			if (used[x][y] == -1) {
				continue;
			} else {
				coordsBW.add(new int[] { x, y });
				used[x][y] = -1;

				revX = xsize - (x + 1);
				revY = ysize - (y + 1);
				x--;
				y--;

				while (revX > 0 && revY > 0 && revM[revX][revY] != 0) {
					actScore = revM[revX][revY];

					if (actScore == revM[revX - 1][revY - 1]
							+ tempScore[xsize - revX][ysize - revY]) {
						if (used[xsize - revX][ysize - revY] != -1) {
							coordsBW.add(new int[] { xsize - revX, ysize - revY });
							used[xsize - revX][ysize - revY] = -1;
							revY--;
							revX--;
						} else {
							break;
						}

					} else if (actScore == revD[revX][revY]) {
						while (revY > 0
								&& revD[revX][revY] == revD[revX][revY - 1]
										+ gapExtend) {
							if (used[xsize - revX][ysize - revY] != -1) {
								coordsBW.add(new int[] { xsize - revX,
										ysize - revY });
								used[xsize - revX][ysize - revY] = -1;
								revY--;
							} else {
								break;
							}
						}
						if (used[xsize - revX][ysize - revY] != -1) {
							coordsBW.add(new int[] { xsize - revX, ysize - revY });
							used[xsize - revX][ysize - revY] = -1;
							revY--;
						} else {
							break;
						}

					} else if (actScore == revI[revX][revY]) {
						while (revX > 0
								&& revI[revX][revY] == revI[revX - 1][revY]
										+ gapExtend) {
							if (used[xsize - revX][ysize - revY] != -1) {
								coordsBW.add(new int[] { xsize - revX,
										ysize - revY });
								used[xsize - revX][ysize - revY] = -1;
								revX--;
							} else {
								break;
							}
						}
						if (used[xsize - revX][ysize - revY] != -1) {
							coordsBW.add(new int[] { xsize - revX, ysize - revY });
							used[xsize - revX][ysize - revY] = -1;
							revX--;
						} else {
							break;
						}
					}
				}

				// do normal traceback
				while (x > 0 && y > 0 && M[x][y] != 0) {
					actScore = M[x][y];

					if (actScore == M[x - 1][y - 1] + tempScore[x][y]) {
						if (used[x][y] != -1) {
							coordsFW.add(new int[] { x, y });
							used[x][y] = -1;
							y--;
							x--;
						} else {
							break;
						}

					} else if (actScore == D[x][y]) {
						while (D[x][y] == D[x][y - 1] + gapExtend && y > 0) {
							if (used[x][y] != -1) {
								coordsFW.add(new int[] { x, y });
								used[x][y] = -1;
								y--;
							} else {
								break;
							}
						}
						if (used[x][y] != -1) {
							coordsFW.add(new int[] { x, y });
							used[x][y] = -1;
							y--;
						} else {
							break;
						}

					} else if (actScore == I[x][y]) {
						while (I[x][y] == I[x - 1][y] + gapExtend && x > 0) {
							if (used[x][y] != -1) {
								coordsFW.add(new int[] { x, y });
								used[x][y] = -1;
								x--;
							} else {
								break;
							}
						}
						if (used[x][y] != -1) {
							coordsFW.add(new int[] { x, y });
							used[x][y] = -1;
							x--;
						} else {
							break;
						}
					}
				}
				Locals mergedali = mergeAlignment(coordsBW, coordsFW);
				if (mergedali.getScore() > maxscore) {
					maxscore = mergedali.getScore();
				}
				result.add(mergedali);
			}
		}
		System.out.print(maxscore + "\n");
		return result;
	}

	private class sortScore implements Comparator<LocalMatch> {

		@Override
		public int compare(LocalMatch arg0, LocalMatch arg1) {
			return arg1.getScore() - arg0.getScore();
		}

	}

	protected class Locals {
		private double evalue;
		private double score;
		private List<int[]> coords = new ArrayList<int[]>();

		public Locals(double score, double evalue, List<int[]> coords) {
			this.evalue = evalue;
			this.coords = coords;
			this.score = score;
		}

		public List<int[]> getCoords() {
			return this.coords;
		}

		public void setEvalue(double evalue) {
			this.evalue = evalue;
		}

		public double getEvalue() {
			return this.evalue;
		}

		public double getScore() {
			return this.score;
		}
	}

	public void closeOut() throws IOException {
		out.close();
	}

	private class LocalMatch {
		private int score;
		private int[] coord = new int[2];

		public LocalMatch(int score, int[] coord) {
			this.score = score;
			this.coord = coord;
		}

		public int[] getCoords() {
			return this.coord;
		}

		public int getScore() {
			return this.score;
		}
	}
}
