package bioinfo.superpos;

import cern.colt.matrix.DoubleFactory1D;
import cern.colt.matrix.DoubleFactory2D;
import cern.colt.matrix.DoubleMatrix1D;
import cern.colt.matrix.DoubleMatrix2D;
import cern.colt.matrix.linalg.Algebra;
import cern.colt.matrix.linalg.SingularValueDecomposition;
import cern.jet.math.Functions;

/**
 * Kabsch class containing calculation routine
 * 
 * @author gobi4
 * 
 */
public class Kabsch {

	/**
	 * 
	 * @param coordinates
	 *            double[][][] containing double[][] of same length
	 * @return Tranfromation calculated by kabsch from P to Q
	 * 
	 *         transformation of Q can be performed by Transformation Objekt
	 *         itself
	 */
	public static Transformation calculateTransformation(
			double[][][] coordinates) {
		// initialize all known quantities
		DoubleFactory2D factory2D = DoubleFactory2D.dense;
		DoubleFactory1D factory1D = DoubleFactory1D.dense;
		DoubleMatrix2D p = factory2D.make(coordinates[0]);
		DoubleMatrix2D q = factory2D.make(coordinates[1]);

		DoubleMatrix1D centroidP = calculateCentroid(p);
		DoubleMatrix1D centroidQ = calculateCentroid(q);

		p = translate(centroidP, p);
		q = translate(centroidQ, q);

		double initError = calculateInitError(p, q);

		Algebra algebra = new Algebra();

		DoubleMatrix2D tP = p.viewDice(); // transpose
		DoubleMatrix2D covarM = algebra.mult(tP, q);

		SingularValueDecomposition svd = new SingularValueDecomposition(covarM);
		DoubleMatrix2D V = svd.getV();
		DoubleMatrix2D S = svd.getS();
		DoubleMatrix2D U = svd.getU();

		// check for right-handed coordinate system
		int d = -1;
		if (algebra.det(V) * algebra.det(U) > 0) {
			d = 1;
		}

		double[][] updateArr = { { 1, 0, 0 }, { 0, 1, 0 }, { 0, 0, d } };
		DoubleMatrix2D updateMat = p.like(3, 3);
		updateMat.assign(updateArr);

		S = algebra.mult(S, updateMat);
		U = algebra.mult(U, updateMat);
		DoubleMatrix2D tU = U.viewDice();

		double err = calcErr(S);

		DoubleMatrix2D rotation = algebra.mult(V, tU);
		double rmsd1 = calcDefRMSD(p, q);
		// no point in calculating the definition of the rmsd since no rotation
		// has been performed
		double rmsd2 = calcEasyRMSD(err, initError, p);
		// if (rmsd1 == rmsd2) {
		// rmsd1 = rmsd2;
		// } else {
		// rmsd1 = -1;
		// }

		DoubleMatrix1D translation = calculateTranslation(centroidQ, centroidP,
				rotation);

		Transformation result = new Transformation(centroidP, centroidQ,
				translation, rotation, rmsd2, -1, -1);
		return result;
	}

	/**
	 * calculates the centroid of a vector, by calculating the mean of every
	 * component of the vector
	 * 
	 * @param vector
	 *            the vector whose centroid we want to compute
	 * @return a new vector containing the centroid coordinates
	 */
	private static DoubleMatrix1D calculateCentroid(DoubleMatrix2D vector) {
		double cx = 0, cy = 0, cz = 0;
		for (int i = 0; i < vector.rows(); i++) {
			cx += vector.get(i, 0);
			cy += vector.get(i, 1);
			cz += vector.get(i, 2);
		}
		cx = cx / (vector.rows() * 1.0);
		cy = cy / (vector.rows() * 1.0);
		cz = cz / (vector.rows() * 1.0);
		double[] result = { cx, cy, cz };
		DoubleFactory1D factory1D = DoubleFactory1D.dense;
		return factory1D.make(result);
	}

	/**
	 * this function translates the coordinates of proteins P and Q by
	 * substituting the centroid
	 * 
	 * @param pCentroid
	 *            the coordinates of protein P's centroid
	 * @param p
	 *            the coordinates of protein P
	 * @param qCentroid
	 *            the coordinates of protein Q's centroid
	 * @param q
	 *            the coordinates of protein Q
	 */
	private static DoubleMatrix2D translate(DoubleMatrix1D pCentroid,
			DoubleMatrix2D p) {
		for (int i = 0; i < p.rows(); i++) {
			p.viewRow(i).assign(pCentroid, Functions.minus);
		}
		return p;
	}

	/**
	 * calculates the initial error
	 * 
	 * @param p
	 *            the coordinates of protein P
	 * @param q
	 *            the coordinates of protein Q
	 * @return
	 */
	private static double calculateInitError(DoubleMatrix2D p, DoubleMatrix2D q) {
		double pInit = 0;
		double qInit = 0;
		double tx = 0, ty = 0, tz = 0;

		for (int i = 0; i < p.rows(); i++) {
			tx = Math.pow(p.get(i, 0), 2);
			ty = Math.pow(p.get(i, 1), 2);
			tz = Math.pow(p.get(i, 2), 2);
			pInit += tx + ty + tz;
		}

		for (int i = 0; i < q.rows(); i++) {
			tx = Math.pow(q.get(i, 0), 2);
			ty = Math.pow(q.get(i, 1), 2);
			tz = Math.pow(q.get(i, 2), 2);
			qInit += tx + ty + tz;
		}
		return pInit + qInit;
	}

	/**
	 * cross-sum of a 3x3 matrix' diagonal
	 * 
	 * @param S
	 *            the matrix in question
	 * @return the sum
	 */
	public static double calcErr(DoubleMatrix2D S) {
		double err = S.get(0, 0) + S.get(1, 1) + S.get(2, 2);
		return err;
	}

	/**
	 * calculate the RMSD of two point groups by the definition
	 * 
	 * @param p
	 *            the coordinates of protein P
	 * @param q
	 *            the coordinates of protein P
	 * @return the RMSD by the definition
	 */
	private static double calcDefRMSD(DoubleMatrix2D p, DoubleMatrix2D q) {
		double rmsd = 0;
		double tx = 0, ty = 0, tz = 0;
		for (int i = 0; i < p.rows(); i++) {
			tx = p.get(i, 0) - q.get(i, 0);
			ty = p.get(i, 1) - q.get(i, 1);
			tz = p.get(i, 2) - q.get(i, 2);

			tx *= tx;
			ty *= ty;
			tz *= tz;

			rmsd = (tx + ty + tz) * (1.0 * p.rows());
		}
		rmsd = Math.sqrt(rmsd);
		return rmsd;
	}

	/**
	 * calculates the RMSD of two point groups via a shortcut
	 * 
	 * @param err
	 *            the sum over the S matrix diagonal
	 * @param e0
	 *            the initial error
	 * @param p
	 *            the coordinates of protein P
	 * @return the RMSD
	 */
	public static double calcEasyRMSD(double err, double e0, DoubleMatrix2D p) {
		double RMSD = Math.sqrt((e0 - 2 * err) / (p.rows() * 1.0));
		return RMSD;
	}

	/**
	 * this function calculates the translation vector needed for the
	 * superposition.
	 * 
	 * @param centroidQ
	 *            the centroid of protein Q's coordinates
	 * @param centroidP
	 *            the centroid of protein P's coordinates
	 * @param rotation
	 *            the rotation matrix
	 * @return the translation vector
	 */
	private static DoubleMatrix1D calculateTranslation(
			DoubleMatrix1D centroidQ, DoubleMatrix1D centroidP,
			DoubleMatrix2D rotation) {
		DoubleMatrix1D translation = centroidP.like();
		centroidQ.assign(Functions.neg);
		// not quite sure that this matrix multiplication should be like this
		multiply1Dwith2D(centroidQ, rotation, translation);

		translation.assign(centroidP, Functions.plus);
		return translation;
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
	private static void multiply1Dwith2D(DoubleMatrix1D mat1D,
			DoubleMatrix2D mat2D, DoubleMatrix1D result) {
		double[] t = new double[3];
		for (int i = 0; i < 3; i++) {
			for (int j = 0; j < 3; j++) {
				t[i] += mat1D.get(j) * mat2D.get(j, i);
			}
		}
		result.assign(t);
	}
}
