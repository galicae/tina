package bioinfo.proteins.fr4gment;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.util.LinkedList;

import bioinfo.proteins.fragm3nt.ProteinFragment;
import bioinfo.proteins.structure.cbeta.CBetaVectorCalculator;
import bioinfo.superpos.Kabsch;
import bioinfo.superpos.Transformation;

/**
 * this class handles all rotations and corresponding garbage for the baseline
 * loop method
 * 
 * @author galicae
 * 
 */
public class LoopRotator {
	private String seq = "ETVQVNLPVSLEDLFVGKKKSFKIGRKGPHGASEKTQIDIQLKPGWKAGTKITYKNQGDYNPQTGRRKTLQFVIQEKS";
	private ProteinFragment template;
	private double[][] querCoord;
	private LinkedList<ProteinFragment> loops;
	private int[] tempPos;
	private int[] curPos;

	public LoopRotator(ProteinFragment template, double[][] querCoord,
			LinkedList<ProteinFragment> loops, int[] tempPos, int[] curPos) {
		super();
		this.template = template;
		this.querCoord = querCoord;
		this.loops = loops;
		this.tempPos = tempPos;
		this.curPos = curPos;
	}

	private ProteinFragment rotateSingleLoop(ProteinFragment curLoop) {
		// first fit loop between gaps in query structure
		// this means superpose the endpoints optimally
		Transformation t = endpointKabsch(curLoop);
		curLoop.setCoordinates(t.transform(curLoop.getAllResidues()));
		double[] referencePoint = curLoop.getResidue(0);

		// now the actual rotation stuff: use huberste's
		// CBetaVectorCalculator(TM)
		// first calculate the vector around which to rotate
		double[] cross = calcRotVector(curLoop);
		LinkedList<ProteinFragment> noClashes = new LinkedList<ProteinFragment>();
		for (int i = 0; i < 60; i++) {
			double[][] trn = CBetaVectorCalculator.calcRotationMatrix(cross,
					5 * i);
			double[][] coordBy5 = curLoop.getAllResidues();
			coordBy5 = CBetaVectorCalculator.matrixMultiplication(
					curLoop.getAllResidues(), trn);
			correctRotation(coordBy5, referencePoint);

			for (int j = 0; j < coordBy5.length; j++) {
				querCoord[curPos[0] + j] = coordBy5[j];
			}

			ProteinFragment tmpResult = new ProteinFragment("test" + i, seq,
					querCoord, querCoord.length);
			if (stericHindrance(curPos, querCoord)) {
				noClashes.add(tmpResult);
			}
		}

		// if (noClashes.size() > 0) {
		// try {
		// BufferedWriter wr = new BufferedWriter(new FileWriter(
		// "rotations" + curLoop.getID()));
		// for (int i = 0; i < noClashes.size(); i++) {
		// wr.write("MODEL        " + (i+1) + "\n");
		// wr.write(noClashes.get(i).toString());
		// wr.write("ENDMDL\n");
		// }
		// wr.close();
		// } catch (Exception e) {
		//
		// }
		// }

		if (noClashes.size() > 0)
			return noClashes.get(0);
		return null;
	}

	public double[][] copyDoubleArray(double[][] a) {
		double[][] b = new double[a.length][a[0].length];
		for (int i = 0; i < a.length; i++) {
			for (int j = 0; j < a[0].length; j++) {
				b[i][j] = a[i][j];
			}
		}
		return b;
	}

	/**
	 * moves rotated matrix so that it coincides with reference point
	 * 
	 * @param rotated
	 * @param reference
	 */
	public void correctRotation(double[][] rotated, double[] reference) {
		double[] correct = new double[reference.length];
		for (int i = 0; i < correct.length; i++) {
			correct[i] = reference[i] - rotated[0][i];
		}

		for (int i = 0; i < rotated.length; i++) {
			for (int j = 0; j < correct.length; j++) {
				rotated[i][j] += correct[j];
			}
		}
	}

	public ProteinFragment rotateAllLoops() {
		LinkedList<ProteinFragment> models = new LinkedList<ProteinFragment>();
		if (loops != null) {
			for (ProteinFragment loop : loops) {
				models.add(rotateSingleLoop(loop));
			}

			if (models.size() > 0)
				return models.get(0);
		}
		
		return null;
	}

	/**
	 * math stuff. Calculation of the rotation vector (normalized length vector
	 * from start to end of current loop fragment)
	 * 
	 * @param curLoop
	 * @return
	 */
	private double[] calcRotVector(ProteinFragment curLoop) {
		double[] cross = new double[3];
		double[] a = curLoop.getResidue(curLoop.getAllResidues().length - 1);
		double[] b = curLoop.getResidue(0);
		cross[0] = a[0] - b[0];
		cross[1] = a[1] - b[1];
		cross[2] = a[2] - b[2];

		double[] origin = { 0, 0, 0 };

		double dist = euclideanDistance(cross, origin);

		cross[0] /= dist;
		cross[1] /= dist;
		cross[2] /= dist;
		return cross;
	}

	private double euclideanDistance(double[] x, double[] y) {
		double result = 0;
		for (int i = 0; i < x.length; i++) {
			result += Math.pow(x[i] - y[i], 2);
		}
		return Math.sqrt(result);
	}

	/**
	 * function to summarize the dull steps for calculating a Kabsch
	 * superposition of the endpoints
	 * 
	 * @param curLoop
	 *            the loop fragment we want to superpose to the ends of the core
	 *            regions surrounding it
	 * @return the {@link Transformation} object corresponding to the move
	 */
	private Transformation endpointKabsch(ProteinFragment curLoop) {
		double[][][] kabschFood = new double[2][2][3];
		kabschFood[1][0] = curLoop.getResidue(0);
		kabschFood[1][1] = curLoop
				.getResidue(curLoop.getAllResidues().length - 1);
		kabschFood[0][0] = template.getResidue(tempPos[0]);
		kabschFood[0][1] = template.getResidue(tempPos[1]);
		return Kabsch.calculateTransformation(kabschFood);
	}

	/**
	 * this function checks a newly inserted loop (its start and end in {@value
	 * position}) for steric clashes
	 * 
	 * @param position
	 *            the start and end of the loop
	 * @param coord
	 *            the coordinates where the loop is inserted
	 * @return true if there are no steric hindrances; false else
	 */
	public boolean stericHindrance(int[] position, double[][] coord) {
		for (int i = position[0]; i < position[1]; i++) {
			for (int j = 0; j < i - 6; j++) {
				if (euclideanDistance(coord[i], coord[j]) < 7)
					return false;
			}
			for (int j = i + 6; j < coord.length; j++) {
				if (euclideanDistance(coord[i], coord[j]) < 7)
					return false;
			}
		}
		return true;
	}

}
