package bioinfo.superpos;

import bioinfo.proteins.PDBEntry;
import cern.colt.matrix.DoubleFactory1D;
import cern.colt.matrix.DoubleMatrix1D;
import cern.colt.matrix.DoubleMatrix2D;
import cern.colt.matrix.impl.DenseDoubleMatrix1D;
import cern.colt.matrix.impl.DenseDoubleMatrix2D;
import cern.jet.math.Functions;

/**
 * value-object class for transformations which are the output of kabsch
 * 
 * @author gobi4
 * 
 */

public class Transformation {

	private DoubleMatrix1D centroidP;
	private DoubleMatrix1D centroidQ;
	private DoubleMatrix1D translation;
	private DoubleMatrix2D rotation;
	private double rmsd;
	private double gdt;
	private double tmscore;

	public Transformation(DoubleMatrix1D centroidP, DoubleMatrix1D centroidQ,
			DoubleMatrix1D translation, DoubleMatrix2D rotation, double rmsd,
			double gdt, double tmscore) {
		super();
		this.centroidP = centroidP;
		this.centroidQ = centroidQ;
		this.translation = translation;
		this.rotation = rotation;
		this.rmsd = rmsd;
		this.gdt = gdt;
		this.tmscore = tmscore;
	}

	/**
	 * transforms Q so that it matches as good as possible on a previously
	 * entered P
	 * 
	 * @param Q
	 *            PDBEntry which was used for kabsch calculation an has to be
	 *            transformed on P now
	 * @return new PDBEntry Q which is transformed old Q
	 */
	public PDBEntry transform(PDBEntry Q) {
		calculateTranslation();
		DoubleFactory1D factory = DoubleFactory1D.dense;
		DoubleMatrix1D tempRow = factory.make(3);

		for (int i = 0; i < Q.length(); i++) {
			for (int j = 0; j < Q.getAminoAcid(i).getAtomNumber(); j++) {
				tempRow = factory.make(Q.getAminoAcid(i).getAtom(j)
						.getPosition());
				multiply1Dwith2D(tempRow, rotation, tempRow);
				tempRow.assign(translation, Functions.plus);
				Q.getAminoAcid(i).getAtom(j).setPosition(tempRow.toArray());
			}
		}
		return Q;
	}

	/**
	 * transforms Q so that it matches as good as possible on a previously
	 * entered P
	 * 
	 * @param Q
	 *            double[][] which was used for kabsch calculation an has to be
	 *            transformed on P now
	 * @return new double[][] Q which is transformed old Q
	 */
	public double[][] transform(double[][] Q) {
		calculateTranslation();
		DoubleFactory1D factory = DoubleFactory1D.dense;
		DoubleMatrix1D tempRow = factory.make(3);

		for (int i = 0; i < Q.length; i++) {
			tempRow = factory.make(Q[i]);
			multiply1Dwith2D(tempRow, rotation, tempRow);
			tempRow.assign(translation, Functions.plus);
			Q[i] = tempRow.toArray();
		}
		return Q;
	}

	/**
	 * this function checks whether the translation was already calculated. If
	 * not, it calculates it with the centroids and the rotation matrix
	 */
	private void calculateTranslation() {
		if (translation != null)
			return;

		centroidQ.assign(Functions.neg);
		// not quite sure that this matrix multiplication should be like this
		multiply1Dwith2D(centroidQ, rotation, translation);

		translation.assign(centroidP, Functions.plus);
	}

	/**
	 * helper function to calculate the cartesian product of a 1D with a 2D
	 * vector
	 * 
	 * @param mat1D
	 *            self explanatory, i think
	 * @param mat2D
	 *            pretty much obvious
	 * @param result
	 *            really?
	 */
	private void multiply1Dwith2D(DoubleMatrix1D mat1D, DoubleMatrix2D mat2D,
			DoubleMatrix1D result) {
		double[] t = new double[3];
		for (int i = 0; i < 3; i++) {
			for (int j = 0; j < 3; j++) {
				t[i] += mat1D.get(j) * mat2D.get(j, i);
			}
		}
		result.assign(t);
	}

	public double getRmsd() {
		return rmsd;
	}

	public double getGdt() {
		return gdt;
	}

	public double getTmscore() {
		return tmscore;
	}
	
	public String getTranslation() {
		return translation.toString();
	}
	
	public String getRotation() {
		return rotation.toString();
	}
	
	public DenseDoubleMatrix1D peekTranslation(){
		return (DenseDoubleMatrix1D)translation;
	}
	
	public DenseDoubleMatrix2D peekRotation(){
		return (DenseDoubleMatrix2D)rotation;
	}
	
	
}
