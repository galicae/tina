package bioinfo.alignment.gotoh;

import java.io.IOException;
import java.io.Writer;
import java.util.Locale;

import bioinfo.Sequence;
import bioinfo.alignment.Alignable;
import bioinfo.alignment.Aligner;
import bioinfo.alignment.Alignment;

/**
 * 
 * @author gobi4 abstract class Gotoh defines interface for alignments if no
 *         method is overridden in extending class it performs an gotoh
 *         alignment on identity matrix with global matrix preparation and
 *         traceback
 * 
 *         for other alignment types the following methods have to be
 *         overridden: + check + align
 */

public abstract class Gotoh implements Aligner {

	public static final int FACTOR = 100;
	protected int gapOpen = 0;
	protected int gapExtend = 0;
	protected int[][] M;
	protected int[][] I;
	protected int[][] D;
	protected Alignable sequence1;
	protected Alignable sequence2;

	/**
	 * Constructor initializing Gotoh with given gap-penalties
	 * 
	 * @param gap
	 *            integer defining gap-costs
	 */
	public Gotoh(double gapOpen, double gapExtend) {
		this.gapOpen = (int) (FACTOR * gapOpen);
		this.gapExtend = (int) (FACTOR * gapExtend);
	}

	@Override
	public Alignment align(Alignable sequence1, Alignable sequence2) {
		this.M = new int[sequence1.length()][sequence2.length()];
		this.I = new int[sequence1.length()][sequence2.length()];
		this.D = new int[sequence1.length()][sequence2.length()];
		this.sequence1 = sequence1;
		this.sequence2 = sequence2;
		return null;
	}

	/**
	 * Override this method in extensions!
	 */
	@Override
	public boolean check(Alignment alignment) {
		return false;
	}

	/**
	 * streams all matrices as tab-separated content
	 * 
	 * @param out
	 *            Writer for example BufferedWriter to file or to System.out
	 */
	public void streamMatricesAsTxt(Writer out) throws IOException {
		out.append("--- score table ---\n\t\t");
		for (int i = 0; i != sequence1.length(); i++) {
			out.append("  " + sequence1.getComp(i) + "\t");
		}
		out.append("\n");
		for (int i = 0; i != sequence2.length(); i++) {
			if (i > 0) {
				out.append(sequence2.getComp(i) + "\t");
			} else {
				out.append(" \t");
			}
			for (int j = 0; j != sequence1.length(); j++) {
				out.append(M[j][i] + "\t");
			}
			out.append("\n");
		}
		out.append("\n");
		out.append("--- deletion table ---\n\t\t");
		for (int i = 0; i != sequence1.length(); i++) {
			out.append("  " + sequence1.getComp(i) + "\t");
		}
		out.append("\n");
		for (int i = 0; i != sequence2.length(); i++) {
			if (i > 0) {
				out.append(sequence2.getComp(i) + "\t");
			} else {
				out.append(" \t");
			}
			for (int j = 0; j != sequence1.length(); j++) {
				out.append(D[j][i] + "\t");
			}
			out.append("\n");
		}
		out.append("\n");
		out.append("--- insertion table ---\n\t");
		for (int i = 0; i != sequence1.length(); i++) {
			out.append("  " + sequence1.getComp(i) + "\t");
		}
		out.append("\n");
		for (int i = 0; i != sequence2.length(); i++) {
			if (i > 0) {
				out.append(sequence2.getComp(i) + "\t");
			} else {
				out.append(" \t");
			}
			for (int j = 0; j != sequence1.length(); j++) {
				out.append(I[j][i] + "\t");
			}
			out.append("\n");
		}
		out.append("\n");
		out.flush();
	}

	/**
	 * streams all matrices as html-formatted content
	 * 
	 * @param out
	 *            Writer, for example BufferedWriter to file or to System.out
	 * 
	 */
	public void streamMatricesAsHtml(Writer out, int factor) throws IOException {
		out.append("<table id=\"score\" style=\"text-align:center;\">\n");
		out.append("<tr><td></td><td></td>");
		for (int i = 0; i != sequence1.length(); i++) {
			out.append("<td><b>" + sequence1.getComp(i) + "</b></td>");
		}
		out.append("</tr>\n");
		for (int i = 0; i <= sequence2.length(); i++) {
			if (i > 0) {
				out.append("<tr><td><b>" + sequence2.getComp(i - 1)
						+ "</b></td>");
			} else {
				out.append("<tr><td></td>");
			}
			for (int j = 0; j <= sequence1.length(); j++) {
				out.append("<td>"
						+ String.format("%6.3f", 1.0d * M[j][i] / factor)
						+ "</td>");
			}
			out.append("</tr>\n");
		}
		out.append("</table>\n");

		out.append("<table id=\"deletion\" style=\"text-align:center;\">\n");
		out.append("<tr><td></td><td></td>");
		for (int i = 0; i != sequence1.length(); i++) {
			out.append("<td><b>" + sequence1.getComp(i) + "</b></td>");
		}
		out.append("</tr>\n");
		for (int i = 0; i <= sequence2.length(); i++) {
			if (i > 0) {
				out.append("<tr><td><b>" + sequence2.getComp(i - 1)
						+ "</b></td>");
			} else {
				out.append("<tr><td></td>");
			}
			for (int j = 0; j <= sequence1.length(); j++) {
				out.append("<td>"
						+ String.format("%6.3f", 1.0d * D[j][i] / factor)
						+ "</td>");
			}
			out.append("</tr>\n");
		}
		out.append("</table>\n");

		out.append("<table id=\"vertical\" style=\"text-align:center;\">\n");
		out.append("<tr><td></td><td></td>");
		for (int i = 0; i != sequence1.length(); i++) {
			out.append("<td><b>" + sequence1.getComp(i) + "</b></td>");
		}
		out.append("</tr>\n");
		for (int i = 0; i <= sequence2.length(); i++) {
			if (i > 0) {
				out.append("<tr><td><b>" + sequence2.getComp(i - 1)
						+ "</b></td>");
			} else {
				out.append("<tr><td></td>");
			}
			for (int j = 0; j <= sequence1.length(); j++) {
				out.append("<td>"
						+ String.format("%6.3f", 1.0d * I[j][i] / factor)
						+ "</td>");
			}
			out.append("</tr>\n");
		}
		out.append("</table>\n");
		out.flush();

	}

	// return M matrix
	public int[][] getM() {
		return this.M;
	}

}
