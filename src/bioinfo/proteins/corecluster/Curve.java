package bioinfo.proteins.corecluster;

import java.io.Serializable;

import cern.colt.matrix.DoubleMatrix1D;
import cern.colt.matrix.linalg.Algebra;
import cern.jet.math.Functions;

public class Curve implements Serializable, Comparable<Curve> {

	/**
	 * 
	 */
	private static final long serialVersionUID = 4584718574380697050L;
	private CompPair<Integer> position;
	private Pair<DoubleMatrix1D> curve;
	private char type;
	private String id;
	
	/**
	 * 
	 * @return sequence length of curve
	 */
	public int getSeqLength() {
		int length = position.getY() - position.getX();
		return length;
	}

	/**
	 * 
	 * @return the length ofthe vector
	 */
	public double getSpaceLength() {
		DoubleMatrix1D x = curve.getX();
		DoubleMatrix1D y = curve.getY();
		x.assign(y, Functions.minus);
		Algebra al = new Algebra();
		double length = Math.sqrt(al.norm2(x));
		
		return length;
	}
	public Curve(String id, CompPair<Integer> position,
			Pair<DoubleMatrix1D> curve, char type) {
		this.id = id;
		this.curve = curve;
		this.position = position;
		this.type = type;
	}

	public CompPair<Integer> getPosition() {
		return position;
	}

	public Pair<DoubleMatrix1D> getCurve() {
		return curve;
	}

	public char getType() {
		return type;
	}

	public String getId() {
		return id;
	}

	@Override
	public int compareTo(Curve o) {
		if (position.getX() > o.getPosition().getX()) {
			return 1;
		} else if (position.getX() < o.getPosition().getX()) {
			return -1;
		} else {
			return 0;
		}
	}

	@Override
	public String toString() {
		return this.position.toString();
	}
}
