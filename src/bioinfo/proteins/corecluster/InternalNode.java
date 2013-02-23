package bioinfo.proteins.corecluster;

import java.util.ArrayList;
import java.util.List;


public class InternalNode {

	private Curve curve;
	private double distance;
	private int length = -1;
	private List<InternalNode> prefixes;
	private InternalNode bestPrefix = null;

	public InternalNode(Curve curve, double distance) {
		this.prefixes = new ArrayList<InternalNode>();
		this.curve = curve;
		this.distance = distance;
	}

	public int getStart() {
		return curve.getPosition().getX();
	}

	public int getEnd() {
		return curve.getPosition().getY();
	}

	public void addPrefix(InternalNode n) {
		prefixes.add(n);
	}

	public int getLength() {
		if (length > 0) {
			return length;
		}
		if (prefixes.isEmpty()) {
			return curve.getPosition().getY() - curve.getPosition().getX();
		}
		int maxLength = -1;
		InternalNode maxLengthNode = null;
		for (InternalNode n : prefixes) {
			if (n.getLength() > maxLength) {
				maxLength = n.getLength();
				maxLengthNode = n;
			}
		}
		bestPrefix = maxLengthNode;
		length = maxLength + (curve.getPosition().getY() - curve.getPosition().getX());
		return length;
	}

	public Curve getCurve() {
		return curve;
	}

	public double getDistance() {
		return distance;
	}

	public InternalNode getBestPrefix() {
		return bestPrefix;
	}

	public String toString() {
		return curve.getPosition().getX() + " " + curve.getPosition().getY() + ": " + length + "\n";
	}
}
