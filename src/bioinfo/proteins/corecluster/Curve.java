package bioinfo.proteins.corecluster;

import java.io.Serializable;

import cern.colt.matrix.DoubleMatrix1D;

public class Curve implements Serializable, Comparable<Curve> {

	/**
	 * 
	 */
	private static final long serialVersionUID = 4584718574380697050L;
	private CompPair<Integer> position;
	private Pair<DoubleMatrix1D> curve;
	private char type;
	private String id;

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
