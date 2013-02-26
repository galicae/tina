package bioinfo.proteins.corecluster;

import java.util.ArrayList;
import java.util.List;

import bioinfo.proteins.fragm3nt.ProteinFragment;

import com.sun.org.apache.xpath.internal.functions.Function;

import cern.colt.function.DoubleDoubleFunction;
import cern.colt.matrix.DoubleMatrix1D;
import cern.colt.matrix.linalg.Algebra;
import cern.jet.math.Functions;

public class MultiCurve {
	private List<Curve> curves;
	private List<ProteinFragment> coordinates;
	private String id;
	private double voronoiArea;
	private double[] theta;

	public MultiCurve(List<Curve> curves, String id, double voronoiArea) {
		this.curves = curves;
		this.id = id;
		this.voronoiArea = voronoiArea;
	}

	public MultiCurve(String id, double voronoiArea) {
		this.curves = new ArrayList<Curve>();
		this.coordinates = new ArrayList<ProteinFragment>();
		this.id = id;
		this.voronoiArea = voronoiArea;
	}

	public Curve[] getCurves() {
		Curve[] result = new Curve[curves.size()];
		for (int i = 0; i < curves.size(); i++)
			result[i] = curves.get(i);
		return result;
	}

	public ProteinFragment[] getStructures() {
		ProteinFragment[] result = new ProteinFragment[coordinates.size()];
		for (int i = 0; i < coordinates.size(); i++)
			result[i] = coordinates.get(i).clone();
		return result;
	}

	public double[][] getCoordinates(int index) {
		try {
			return coordinates.get(index).getAllResidues();
		} catch (NullPointerException e) {
			return null;
		}
	}

	public void addElement(Curve c, ProteinFragment f) {
		curves.add(c);
		coordinates.add(f);
	}

	public double getVoronoiArea() {
		return voronoiArea;
	}

	public String getId() {
		return id;
	}

	/**
	 * this function calculates the angle between the curves in the indices
	 * 
	 * @param curve1Index
	 *            the index of the first curve (starting at 0)
	 * @param curve2Index
	 *            the index of the second curve (starting at 0)
	 * @return
	 */
	public double calculateAngle(int curve1Index, int curve2Index) {
		double theta = 0;
		Curve curve1 = curves.get(curve1Index);
		Curve curve2 = curves.get(curve2Index);

		DoubleMatrix1D x1 = curve1.getCurve().getX();
		DoubleMatrix1D y1 = curve1.getCurve().getY();
		DoubleMatrix1D x2 = curve2.getCurve().getX();
		DoubleMatrix1D y2 = curve2.getCurve().getY();

		x1.assign(y1, Functions.minus);
		x2.assign(y2, Functions.minus);

		Algebra al = new Algebra();

		double up = al.mult(x1, x2);

		double a = Math.sqrt(al.norm2(x1));
		double b = Math.sqrt(al.norm2(x2));
		double down = a * b;
		theta = Math.acos(up / down);

		return theta;
	}

	/**
	 * calculates the angles between all curves, and saves them in the theta
	 * array
	 */
	public void calculateAllAngles() {
		theta = new double[((curves.size() - 1) * curves.size()) / 2];
		for (int i = 0; i < curves.size() - 1; i++) {
			for (int j = i + 1; j < curves.size(); j++) {
				theta[((j * (j - 1)) / 2) + i] = calculateAngle(i, j);
			}
		}
	}

	public double getTheta(int i, int j) {
		int index = ((j * (j - 1)) / 2) + i;
		return theta[index];
	}
}
