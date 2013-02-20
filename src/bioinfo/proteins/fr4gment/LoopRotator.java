package bioinfo.proteins.fr4gment;

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
	private ProteinFragment template;
	private double[][] querCoord;
	private LinkedList<ProteinFragment> loops;
	private int[] tempPos;

	public LoopRotator(ProteinFragment template, double[][] querCoord,
			LinkedList<ProteinFragment> loops, int[] tempPos) {
		super();
		this.template = template;
		this.querCoord = querCoord;
		this.loops = loops;
		this.tempPos = tempPos;
	}

	public double rotateSingleLoop(ProteinFragment curLoop) {
		// first fit loop between gaps in query structure
		// this means superpose the endpoints optimally
		Transformation t = endpointKabsch(curLoop);
		curLoop.setCoordinates(t.transform(curLoop.getAllResidues()));

		// now the actual rotation stuff: use huberste's
		// CBetaVectorCalculator(TM)
		CBetaVectorCalculator.calcRotationMatrix(cross, phi);
		
		return -1;
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
		kabschFood[0][0] = curLoop.getResidue(0);
		kabschFood[0][1] = curLoop
				.getResidue(curLoop.getAllResidues().length - 1);
		kabschFood[1][0] = template.getResidue(tempPos[0]);
		kabschFood[1][1] = template.getResidue(tempPos[1]);
		return Kabsch.calculateTransformation(kabschFood);
	}

	public boolean stericHindrance(int[] position) {

		return false;
	}

}
