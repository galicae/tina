/******************************************************************************
 * huberdp.scoring.MagicScoring.java                                          *
 *                                                                            *
 * Contains the MagicScoring class which uses TMAlign as an scoring indicator *
 *                                                                            *
 * This file is best read at line width 80 and tab width 4.                   *
 *                                                                   huberste *
 ******************************************************************************/
package huberdp.scoring;

import static util.Util.*;

import java.util.LinkedList;

import bioinfo.alignment.Threading;
import bioinfo.proteins.AminoAcid;
import bioinfo.proteins.Atom;
import bioinfo.proteins.PDBEntry;
import bioinfo.proteins.PDBFileReader;
import bioinfo.superpos.tmalign.TMAlign;
import huberdp.Scoring;

/**
 * @author huberste
 * @lastchange 2013-02-24
 */
public class MagicScoring implements Scoring {

	/**
	 * Strings for paths
	 */
	String pdbpath;
	
	String tmppath;

	TMAlign tmalign;

	/**
	 * scoring matrix
	 */
	double[][] matrix;

	/**
	 * IDs (xxxxA00) of template and target
	 */
	String template, target;

	/**
	 * constructor for a MagicScoring
	 */
	public MagicScoring(String pdbpath, String tmpath, String tmppath) {
		this.pdbpath = pdbpath;
		this.tmppath = tmppath;
		this.tmalign = new TMAlign(tmpath);
	}

	/**
	 * initializes stuff.
	 * 
	 * @param template
	 * @param target
	 */
	private void init(String template, String target) {
		// execute TMAlign
		this.template = template;
		this.target = target;
		String templatepdbfile = pdbpath + template + ".pdb";
		String targetpdbfile = pdbpath + target + ".pdb";
		// calculate Rotation Matrix
		tmalign.align(templatepdbfile, targetpdbfile, tmppath+"rmt");

		// load rotation matrix file
		double[][] u = TMAlign.loadRotationMatrixFromFile(tmppath+"rmt");

		// rotate the template structure
		PDBEntry templateStruct = new PDBFileReader()
				.readPDBFromFile(templatepdbfile);
		PDBEntry targetStruct = new PDBFileReader()
				.readPDBFromFile(targetpdbfile);
		LinkedList<AminoAcid> aminos = new LinkedList<AminoAcid>();

		// rotate template to target
		// for each AA in template
		for (int i = 0; i < templateStruct.length(); i++) {
			LinkedList<Atom> atoms = new LinkedList<Atom>();
			// for each Atom in AA
			for (int j = 0; j < templateStruct.getComp(i).getAtomNumber(); j++) {
				double[] v = templateStruct.getComp(i).getAtom(j).getPosition();
				double[] w = new double[3];
				w[0] = u[0][0] + u[1][0]*v[0] + u[2][0]*v[1] + u[3][0]*v[2];
				w[1] = u[0][1] + u[1][1]*v[0] + u[2][1]*v[1] + u[3][1]*v[2];
				w[2] = u[0][2] + u[1][2]*v[0] + u[2][2]*v[1] + u[3][2]*v[2];
				atoms.add(new Atom(templateStruct.getComp(i).getAtom(j)
						.getType(), w));
			}
			aminos.add(new AminoAcid(templateStruct.getComp(i).getName(),
					templateStruct.getComp(i).getResIndex(), atoms
							.toArray(new Atom[0])));
		}
		
		// calculate Matrix!
		matrix = new double[templateStruct.length()][targetStruct.length()];

		for (int m = 0; m < templateStruct.length(); m++) {
			for (int n = 0; n < targetStruct.length(); n++) {
				matrix[m][n] = calcDistance(targetStruct.getComp(n),
						aminos.get(m));
			}
		}
	}

	@Override
	public double score(Threading threading) {
		if (target == null || template == null
				|| (!target.equals(threading.getSequence().getId()))
				|| (!template.equals(threading.getStructure().getID()))) {
			init(threading.getStructure().getID(), threading.getSequence()
					.getId());
		}

		double result = 0.0;
		int temppos = 0;
		int targpos = 0;
		int[][] rows = threading.getRows();
		for (int i = 0; i < rows[0].length; i++) {
			if (rows[0][i] == -1) { // insertion
				result += getInsertionScore(threading, targpos);
				targpos++;
			} else if (rows[1][i] == -1) { // deletion
				result += getDeletionScore(threading, temppos);
				temppos++;
			} else { // match
				result += getScore(threading, temppos, targpos);
				temppos++;
				targpos++;
			}
		}
		return result;
	}

	@Override
	public double getScore(Threading t, int m, int n) {
		if (template == null || target == null
				|| (!template.equals(t.getStructure().getID()))
				|| (!target.equals(t.getSequence().getId()))) {
			init(t.getStructure().getID(), t.getSequence().getId());
		}

		if (matrix[m][n] < 1)
			return 10;
		else if (matrix[m][n] < 2)
			return 9;
		else if (matrix[m][n] < 3)
			return 8;
		else if (matrix[m][n] < 4)
			return 7;
		else if (matrix[m][n] < 5)
			return 6;
		else if (matrix[m][n] < 6)
			return 5;
		else if (matrix[m][n] < 7)
			return 3;
		else if (matrix[m][n] < 8)
			return 1;
		else
			return 0;
	}

	@Override
	public double getInsertionScore(Threading t, int n) {
		return -5.0;
	}

	@Override
	public double getDeletionScore(Threading t, int m) {
		return -5.0;
	}

}

/******************************************************************************
 * "A question that sometimes drives me hazy:                                 *
 *  Am I or are the others crazy?"                                            *
 * - Albert Einstein (1879 - 1955)                                            *
 ******************************************************************************/
