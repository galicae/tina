package bioinfo.proteins.corecluster;

import bioinfo.proteins.AminoAcid;
import bioinfo.proteins.AtomType;
import bioinfo.proteins.PDBEntry;

import cern.colt.matrix.DoubleMatrix1D;
import cern.colt.matrix.DoubleMatrix2D;
import cern.colt.matrix.impl.DenseDoubleMatrix1D;
import cern.colt.matrix.impl.DenseDoubleMatrix2D;

public class VectorMath {

	private static final int X = 0, Y = 1, Z = 2;

	public static Curve getCurve(PDBEntry chain, CompPair<Integer> position, char sse)  {
		DoubleMatrix2D elementPoints = getSubElementPointsBB(chain, position);
		DoubleMatrix2D distancePoints = getSubDistancePoints(chain, position);
		Pair<DoubleMatrix1D> line = getLine(distancePoints);
		DoubleMatrix1D start = VectorMath.getPoint(line, VectorMath.getScalar(line, elementPoints.viewRow(0)));
		DoubleMatrix1D end = VectorMath.getPoint(line, VectorMath.getScalar(line, elementPoints.viewRow(elementPoints.rows() - 1)));
		return new Curve(chain.getID(), position, new Pair<DoubleMatrix1D>(start, end), sse);
	}

	public static DoubleMatrix2D getSubElementPointsBB(PDBEntry c, CompPair<Integer> position)  {
		DoubleMatrix2D result = new DenseDoubleMatrix2D((position.getY() - position.getX()) * 4, 3);
		for (int index = position.getX(); index < position.getY(); index++) {
			result.set((index - position.getX()) * 4, X, getCoords(c.getAminoAcid(index), "N").get(X));
			result.set((index - position.getX()) * 4, Y, getCoords(c.getAminoAcid(index), "N").get(Y));
			result.set((index - position.getX()) * 4, Z, getCoords(c.getAminoAcid(index), "N").get(Z));
			result.set((index - position.getX()) * 4 + 1, X, getCoords(c.getAminoAcid(index), "CA").get(X));
			result.set((index - position.getX()) * 4 + 1, Y, getCoords(c.getAminoAcid(index), "CA").get(Y));
			result.set((index - position.getX()) * 4 + 1, Z, getCoords(c.getAminoAcid(index), "CA").get(Z));
			result.set((index - position.getX()) * 4 + 2, X, getCoords(c.getAminoAcid(index), "C").get(X));
			result.set((index - position.getX()) * 4 + 2, Y, getCoords(c.getAminoAcid(index), "C").get(Y));
			result.set((index - position.getX()) * 4 + 2, Z, getCoords(c.getAminoAcid(index), "C").get(Z));
			result.set((index - position.getX()) * 4 + 3, X, getCoords(c.getAminoAcid(index), "O").get(X));
			result.set((index - position.getX()) * 4 + 3, Y, getCoords(c.getAminoAcid(index), "O").get(Y));
			result.set((index - position.getX()) * 4 + 3, Z, getCoords(c.getAminoAcid(index), "O").get(Z));
		}
		return result;
	}

	public static DoubleMatrix2D getSubElementPointsCA(PDBEntry c, CompPair<Integer> position)  {
		DoubleMatrix2D result = new DenseDoubleMatrix2D((position.getY() - position.getX()), 3);
		for (int index = position.getX(); index < position.getY(); index++) {
			result.set(index - position.getX(), X, getCoords(c.getAminoAcid(index), "CA").get(X));
			result.set(index - position.getX(), Y, getCoords(c.getAminoAcid(index), "CA").get(Y));
			result.set(index - position.getX(), Z, getCoords(c.getAminoAcid(index), "CA").get(Z));
		}
		return result;
	}

	/**
	 * returns the vectors of the ca-atom coordinates as a matrix
	 * 
	 * @param c
	 *            protein chain with residue groups of ca-atoms
	 * @param position
	 *            start and end point of the residues
	 * @return matrix with row vectors of ca-atom coords
	 * @
	 */
	public static DoubleMatrix2D getSubDistancePoints(PDBEntry c, CompPair<Integer> position)  {
		DoubleMatrix2D template = getSubElementPointsBB(c, position);
		DoubleMatrix2D target = getDistancePoints(template);
		Calculator calculator = new Calculator(template, target);
		return calculator.getSuperposition();
	}

	/**
	 * returns the coordinates of the ca-atom in the group as a vector
	 * 
	 * @param group
	 *            group with ca-atom
	 * @return vector of ca-atom coords
	 * @
	 */
	private static DoubleMatrix1D getCoords(AminoAcid a, String atomName) {
		return new DenseDoubleMatrix1D(a.getAtomByType(AtomType.createFromString(atomName)).getPosition());
		
	}

	/**
	 * Calculates a set of points on the z-axis based on the distances between
	 * the points in data
	 * 
	 * @param data
	 * @return line on the z-axis
	 */
	public static DoubleMatrix2D getDistancePoints(DoubleMatrix2D data) {
		DoubleMatrix2D result = new DenseDoubleMatrix2D(data.rows(), data.columns());
		double zValue = 0.0;
		for (int row = 1; row < data.rows(); row++) {
			zValue += Calculator.getDistance(data.viewRow(row - 1), data.viewRow(row));
			result.set(row, X, 0.0);
			result.set(row, Y, 0.0);
			result.set(row, Z, zValue);
		}
		return result;
	}

	/**
	 * checks if a curve fits a set of points
	 * 
	 * @param curve
	 *            a pair of points defining a curve
	 * @param points
	 *            a set of points
	 * @param avg
	 *            arithmetic mean
	 * @param sd
	 *            standard deviation
	 * @return true if there is no deviation from avg+-sd
	 * @
	 */
	public static boolean isVectorNormal(Curve curve, PDBEntry chain, double avg, double sd)  {
		return (getAverageDistance(curve, chain) <= (avg + sd) && getAverageDistance(curve, chain) >= (avg - sd));
		// return (getMaxDistance(curve, chain).getY() < (avg + 2 * sd) &&
		// getMinDistance(curve, chain).getY() > (avg - 2 * sd));
	}

	/**
	 * calculates the scalar for a given line so that the distance from the line
	 * to the point is minimal
	 * 
	 * @param line
	 *            a pair of vectors, first support vector, second orientation
	 *            vector
	 * @param point
	 *            another point
	 * @return scaling factor for the line
	 */
	public static double getScalar(Pair<DoubleMatrix1D> line, DoubleMatrix1D point) {
		double mx = line.getY().get(X), my = line.getY().get(Y), mz = line.getY().get(Z);
		double ax = point.get(X), ay = point.get(Y), az = point.get(Z);
		double bx = line.getX().get(X), by = line.getX().get(Y), bz = line.getX().get(Z);
		return (mx * (ax - bx) + my * (ay - by) + mz * (az - bz)) / (mx * mx + my * my + mz * mz);
	}

	/**
	 * calculates a point on a line given a factor for scaling the line
	 * 
	 * @param line
	 *            a pair of vectors, first support vector, second orientation
	 *            vector
	 * @param scalar
	 *            scaling factor for the line
	 * @return point on the line
	 */
	public static DoubleMatrix1D getPoint(Pair<DoubleMatrix1D> line, double scalar) {
		DoubleMatrix1D result = new DenseDoubleMatrix1D(3);
		result.set(X, line.getX().get(X) + line.getY().get(X) * scalar);
		result.set(Y, line.getX().get(Y) + line.getY().get(Y) * scalar);
		result.set(Z, line.getX().get(Z) + line.getY().get(Z) * scalar);
		return result;
	}

	/**
	 * calculates a line crossing a pair of vectors
	 * 
	 * @param curve
	 *            a pair of points
	 * @return pair of vectors, first support vector, second orientation vector
	 */
	public static Pair<DoubleMatrix1D> getLine(Pair<DoubleMatrix1D> curve) {
		DoubleMatrix1D supportVector = new DenseDoubleMatrix1D(3);
		supportVector.set(X, curve.getX().get(X));
		supportVector.set(Y, curve.getX().get(Y));
		supportVector.set(Z, curve.getX().get(Z));
		DoubleMatrix1D orientationVector = new DenseDoubleMatrix1D(3);
		orientationVector = Calculator.add(curve.getY(), Calculator.invertValues(supportVector));
		return new Pair<DoubleMatrix1D>(supportVector, orientationVector);
	}

	/**
	 * calculates the minimal distance of a point to a curve, the curve is first
	 * converted to a line
	 * 
	 * @param curve
	 *            a pair of points defining a curve
	 * @param point
	 *            a point different from the two others
	 * @return minimal distance between the curve and the point
	 */
	public static double getDistancePointToLine(Pair<DoubleMatrix1D> curve, DoubleMatrix1D point) {
		return Calculator.getDistance(point, getPoint(getLine(curve), getScalar(getLine(curve), point)));
	}

	/**
	 * calculates the maximum distance from a set of points and a curve defined
	 * by two points
	 * 
	 * @param curve
	 *            a pair of points defining a curve
	 * @param points
	 *            a set of points
	 * @return a pair of values, first the position of the maximum in the set,
	 *         second the maximum
	 * @
	 */
	public static Pair<Double> getMaxDistance(Curve curve, PDBEntry chain)  {
		double max = Double.MIN_VALUE, pos = 0.0;
		for (int i = curve.getPosition().getX(); i < curve.getPosition().getY(); i++) {
			double dist = getDistancePointToLine(curve.getCurve(), getCoords(chain.getAminoAcid(i), "CA"));
			if (dist > max) {
				max = dist;
				pos = i;
			}
		}
		return new Pair<Double>(pos, max);
	}

	/**
	 * calculates the minimum distance from a set of points and a curve defined
	 * by two points
	 * 
	 * @param curve
	 *            a pair of points defining a curve
	 * @param points
	 *            a set of points
	 * @return a pair of values, first the position of the minimum in the set,
	 *         second the minimum
	 * @
	 */
	public static Pair<Double> getMinDistance(Curve curve, PDBEntry chain)  {
		double min = Double.MAX_VALUE, pos = 0.0;
		for (int i = curve.getPosition().getX(); i < curve.getPosition().getY(); i++) {
			double dist = getDistancePointToLine(curve.getCurve(), getCoords(chain.getAminoAcid(i), "CA"));
			if (dist < min) {
				min = dist;
				pos = i;
			}
		}
		return new Pair<Double>(pos, min);
	}

	/**
	 * calculates the average distance between a curve and a set of points
	 * 
	 * @param curve
	 *            curve defined by two points
	 * @param points
	 *            a set of points
	 * @return average distance
	 * @
	 */
	public static double getAverageDistance(Curve curve, PDBEntry chain)  {
		double result = 0.0;
		int counter = 0;
		for (int i = curve.getPosition().getX(); i < curve.getPosition().getY(); i++) {
			result += getDistancePointToLine(curve.getCurve(), getCoords(chain.getAminoAcid(i), "CA"));
			counter++;
		}
		return result / counter;
	}

	/**
	 * calculates the line between the first and the last points of the set
	 * 
	 * @param points
	 *            a set of points
	 * @return pair of vectors, first the support vector, second the orientation
	 *         vector
	 */
	public static Pair<DoubleMatrix1D> getLine(DoubleMatrix2D points) {
		DoubleMatrix1D supportVector = new DenseDoubleMatrix1D(3);
		supportVector.set(X, points.get(0, X));
		supportVector.set(Y, points.get(0, Y));
		supportVector.set(Z, points.get(0, Z));
		DoubleMatrix1D orientationVector = new DenseDoubleMatrix1D(3);
		orientationVector = Calculator.add(points.viewRow(points.rows() - 1), Calculator.invertValues(supportVector));
		return new Pair<DoubleMatrix1D>(supportVector, orientationVector);
	}
}
