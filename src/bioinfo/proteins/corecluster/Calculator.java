package bioinfo.proteins.corecluster;

import cern.colt.matrix.DoubleMatrix1D;
import cern.colt.matrix.DoubleMatrix2D;
import cern.colt.matrix.impl.DenseDoubleMatrix1D;
import cern.colt.matrix.impl.DenseDoubleMatrix2D;
import cern.colt.matrix.linalg.Algebra;
import cern.colt.matrix.linalg.SingularValueDecomposition;

/**
 * 
 * @author ike
 */
public class Calculator {

	private static final int X = 0, Y = 1, Z = 2;
	private static Algebra aAlgebra = new Algebra();
	private static double[][] aIdentityMatrix = { { 1, 0, 0 }, { 0, 1, 0 }, { 0, 0, 1 } };
	private DoubleMatrix2D aTemplate, aTarget, aRotation = null, aSuperposition = null;
	private DoubleMatrix1D aOrigoVectorTemplate = null, aOrigoVectorTarget = null, aTranslationVector = null;
	private double aKabschRMSD = Double.NaN;

	public Calculator(DoubleMatrix2D pTemplate, DoubleMatrix2D pTarget) {
		aTemplate = pTemplate;
		aTarget = pTarget;
	}

	public DoubleMatrix2D getRotationMatrix() {
		if (aRotation == null) {
			DoubleMatrix2D cP = translateToOrigo(aTemplate);
			DoubleMatrix2D cQ = translateToOrigo(aTarget);
			double initError = getInitError(cP, cQ);
			DoubleMatrix2D tcP = aAlgebra.transpose(cP);
			DoubleMatrix2D covarM = aAlgebra.mult(tcP, cQ);
			SingularValueDecomposition svd = new SingularValueDecomposition(covarM);
			DoubleMatrix2D v = svd.getV();
			DoubleMatrix2D s = svd.getS();
			DoubleMatrix2D u = svd.getU();
			aIdentityMatrix[2][2] = Math.signum(aAlgebra.det(v) * aAlgebra.det(aAlgebra.transpose(u)));
			s = aAlgebra.mult(s, new DenseDoubleMatrix2D(aIdentityMatrix));
			u = aAlgebra.mult(u, new DenseDoubleMatrix2D(aIdentityMatrix));
			aKabschRMSD = Math.sqrt(Math.abs(initError - 2 * getDiagonalSum(s)) / aTemplate.rows());
			aRotation = aAlgebra.mult(v, aAlgebra.transpose(u));
		}
		return aRotation;
	}

	public DoubleMatrix2D getTemplate() {
		return aTemplate;
	}

	public DoubleMatrix2D getTarget() {
		return aTarget;
	}

	public DoubleMatrix1D getTranslationMatrix() {
		if (aTranslationVector == null) {
			aTranslationVector = add(mult(invertValues(getTargetOrigoVector()), getRotationMatrix()), getTemplateOrigoVector());
		}
		return aTranslationVector;
	}

	public DoubleMatrix1D getTemplateOrigoVector() {
		if (aOrigoVectorTemplate == null) {
			aOrigoVectorTemplate = getOrigoTranslationVector(aTemplate);
		}
		return aOrigoVectorTemplate;
	}

	public DoubleMatrix1D getTargetOrigoVector() {
		if (aOrigoVectorTarget == null) {
			aOrigoVectorTarget = getOrigoTranslationVector(aTarget);
		}
		return aOrigoVectorTarget;
	}

	public double getKabschRMSD() {
		if (Double.isNaN(aKabschRMSD)) {
			getRotationMatrix();
		}
		return aKabschRMSD;
	}

	public DoubleMatrix2D getSuperposition() {
		if (aSuperposition == null) {
			DoubleMatrix2D result = new DenseDoubleMatrix2D(aTemplate.rows(), aTemplate.columns());
			DoubleMatrix1D current;
			for (int row = 0; row < aTarget.rows(); row++) {
				current = mult(aTarget.viewRow(row), getRotationMatrix());
				result.set(row, X, current.get(X));
				result.set(row, Y, current.get(Y));
				result.set(row, Z, current.get(Z));
			}
			aSuperposition = add(result, getTranslationMatrix());
		}
		return aSuperposition;
	}

	public DoubleMatrix2D superpositionate(DoubleMatrix2D pFullStructure) {
		DoubleMatrix2D result = new DenseDoubleMatrix2D(pFullStructure.rows(), pFullStructure.columns());
		DoubleMatrix1D current;
		for (int row = 0; row < pFullStructure.rows(); row++) {
			current = mult(pFullStructure.viewRow(row), getRotationMatrix());
			result.set(row, X, current.get(X));
			result.set(row, Y, current.get(Y));
			result.set(row, Z, current.get(Z));
		}
		return add(result, getTranslationMatrix());
	}
	
	public DoubleMatrix1D superpositionate(DoubleMatrix1D pPoint) {
		DoubleMatrix1D result = new DenseDoubleMatrix1D(3);
		result = mult(pPoint, getRotationMatrix());
		return add(result, getTranslationMatrix());
	}

	public static int[][] getDistances(DoubleMatrix2D pFirstMatrix, DoubleMatrix2D pSecondMatrix, int pShift) {
		int[][] result = new int[pFirstMatrix.rows()][pSecondMatrix.rows()];
		for (int firstRow = 0; firstRow < pFirstMatrix.rows(); firstRow++) {
			for (int secondRow = 0; secondRow < pSecondMatrix.rows(); secondRow++) {
				result[firstRow][secondRow] = (int) ((pShift - getDistance(pFirstMatrix.viewRow(firstRow), pSecondMatrix.viewRow(secondRow))) * 10000);
			}
		}
		return result;
	}

	private static DoubleMatrix1D mult(DoubleMatrix1D pVector, DoubleMatrix2D pMatrix) {
		DoubleMatrix1D result = new DenseDoubleMatrix1D(pVector.size());
		double x = pVector.get(X), y = pVector.get(Y), z = pVector.get(Z);
		result.set(X, x * pMatrix.get(X, X) + y * pMatrix.get(Y, X) + z * pMatrix.get(Z, X));
		result.set(Y, x * pMatrix.get(X, Y) + y * pMatrix.get(Y, Y) + z * pMatrix.get(Z, Y));
		result.set(Z, x * pMatrix.get(X, Z) + y * pMatrix.get(Y, Z) + z * pMatrix.get(Z, Z));
		return result;
	}

	public static DoubleMatrix1D add(DoubleMatrix1D pMatrix1, DoubleMatrix1D pMatrix2) {
		DoubleMatrix1D result = new DenseDoubleMatrix1D(pMatrix1.size());
		result.set(X, pMatrix1.get(X) + pMatrix2.get(X));
		result.set(Y, pMatrix1.get(Y) + pMatrix2.get(Y));
		result.set(Z, pMatrix1.get(Z) + pMatrix2.get(Z));
		return result;
	}

	public static DoubleMatrix1D invertValues(DoubleMatrix1D pMatrix1D) {
		DoubleMatrix1D result = new DenseDoubleMatrix1D(pMatrix1D.size());
		result.set(X, -pMatrix1D.get(X));
		result.set(Y, -pMatrix1D.get(Y));
		result.set(Z, -pMatrix1D.get(Z));
		return result;
	}

	private static DoubleMatrix2D add(DoubleMatrix2D pMatrix1, DoubleMatrix1D pMatrix2) {
		for (int row = 0; row < pMatrix1.rows(); row++) {
			pMatrix1.set(row, X, pMatrix1.get(row, X) + pMatrix2.get(X));
			pMatrix1.set(row, Y, pMatrix1.get(row, Y) + pMatrix2.get(Y));
			pMatrix1.set(row, Z, pMatrix1.get(row, Z) + pMatrix2.get(Z));
		}
		return pMatrix1;
	}

	private static DoubleMatrix1D getOrigoTranslationVector(DoubleMatrix2D pMatrix2D) {
		DoubleMatrix1D translationVector = new DenseDoubleMatrix1D(3);
		translationVector.set(X, rowSum(pMatrix2D, X) / pMatrix2D.rows());
		translationVector.set(Y, rowSum(pMatrix2D, Y) / pMatrix2D.rows());
		translationVector.set(Z, rowSum(pMatrix2D, Z) / pMatrix2D.rows());
		return translationVector;
	}

	private static DoubleMatrix2D translateToOrigo(DoubleMatrix2D pMatrix2D) {
		return subtractFromAllRows(pMatrix2D, getOrigoTranslationVector(pMatrix2D));
	}

	private static double getInitError(DoubleMatrix2D pQ, DoubleMatrix2D pP) {
		double result = 0.0;
		for (int row = 0; row < pQ.rows(); row++) {
			result += Math.pow(pQ.get(row, X), 2.0) + Math.pow(pQ.get(row, Y), 2.0) + Math.pow(pQ.get(row, Z), 2.0);
		}
		for (int row = 0; row < pP.rows(); row++) {
			result += Math.pow(pP.get(row, X), 2.0) + Math.pow(pP.get(row, Y), 2.0) + Math.pow(pP.get(row, Z), 2.0);
		}
		return result;
	}

	private static double rowSum(DoubleMatrix2D pMatrix2D, int pColIndex) {
		double result = 0.0;
		for (int row = 0; row < pMatrix2D.rows(); row++) {
			result += pMatrix2D.get(row, pColIndex);
		}
		return result;
	}

	private static DoubleMatrix2D subtractFromAllRows(DoubleMatrix2D pMatrix2D, DoubleMatrix1D pSubtractionVector) {
		DoubleMatrix2D result = new DenseDoubleMatrix2D(pMatrix2D.rows(), pMatrix2D.columns());
		for (int row = 0; row < pMatrix2D.rows(); row++) {
			result.set(row, X, pMatrix2D.get(row, X) - pSubtractionVector.get(X));
			result.set(row, Y, pMatrix2D.get(row, Y) - pSubtractionVector.get(Y));
			result.set(row, Z, pMatrix2D.get(row, Z) - pSubtractionVector.get(Z));
		}
		return result;
	}

	private static double getDiagonalSum(DoubleMatrix2D pMatrix2D) {
		double result = 0.0;
		for (int i = 0; i < pMatrix2D.rows(); i++) {
			result += pMatrix2D.get(i, i);
		}
		return result;
	}

	public static double getRMSD(DoubleMatrix2D pMatrix1, DoubleMatrix2D pMatrix2) {
		double result = 0.0;
		for (int row = 0; row < pMatrix1.rows(); row++) {
			result += Math.pow(pMatrix1.get(row, X) - pMatrix2.get(row, X), 2.0) + Math.pow(pMatrix1.get(row, Y) - pMatrix2.get(row, Y), 2.0)
					+ Math.pow(pMatrix1.get(row, Z) - pMatrix2.get(row, Z), 2.0);
		}
		return Math.sqrt(result / pMatrix1.rows());
	}

	public static double getDistance(DoubleMatrix1D pPoint1, DoubleMatrix1D pPoint2) {
		return Math.sqrt(Math.pow(pPoint1.get(X) - pPoint2.get(X), 2.0) + Math.pow(pPoint1.get(Y) - pPoint2.get(Y), 2.0) + Math.pow(pPoint1.get(Z) - pPoint2.get(Z), 2.0));
	}

	public static double getGDT_TSScore(DoubleMatrix2D pTemplate, DoubleMatrix2D pTarget) {
		int counter1 = 0, counter2 = 0, counter4 = 0, counter8 = 0;
		double distance;
		for (int row = 0; row < pTemplate.rows(); row++) {
			distance = getDistance(pTemplate.viewRow(row), pTarget.viewRow(row));
			if (distance < 1.0) {
				counter1++;
			}
			if (distance < 2.0) {
				counter2++;
			}
			if (distance < 4.0) {
				counter4++;
			}
			if (distance < 8.0) {
				counter8++;
			}
		}
		return (counter1 + counter2 + counter4 + counter8) / (4.0 * pTemplate.rows());
	}

}
