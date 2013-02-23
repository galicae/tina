package bioinfo.proteins.corecluster;

import java.util.ArrayList;
import java.util.List;

import com.sun.org.apache.xpath.internal.functions.Function;

import cern.colt.function.DoubleDoubleFunction;
import cern.colt.matrix.DoubleMatrix1D;
import cern.colt.matrix.linalg.Algebra;
import cern.jet.math.Functions;

public class MultiCurve {
	private List<Curve> curves;
	private String id;
	private double voronoiArea;

	public MultiCurve(List<Curve> curves, String id, double voronoiArea) {
		this.curves = curves;
		this.id = id;
		this.voronoiArea = voronoiArea;
	}
	
	public MultiCurve(String id, double voronoiArea) {
		this.curves = new ArrayList<Curve>();
		this.id = id;
		this.voronoiArea = voronoiArea;
	}

	public Curve[] getCurves() {
		Curve[] result = new Curve[curves.size()];
		for (int i = 0; i < curves.size(); i++)
			result[i] = curves.get(i);
		return result;
	}
	
	public void addElement(Curve c) {
		curves.add(c);
	}

	public double getVoronoiArea() {
		return voronoiArea;
	}

	public String getId() {
		return id;
	}
	
	public double calculateAngle(int curve1Index, int curve2Index) {
		double theta = 0;
		Curve curve1 = curves.get(curve1Index);
		Curve curve2 = curves.get(curve2Index);
		
		DoubleMatrix1D x1 = curve1.getCurve().getX();
		DoubleMatrix1D y1 = curve1.getCurve().getY();
		DoubleMatrix1D x2 = curve2.getCurve().getX();
		DoubleMatrix1D y2 = curve2.getCurve().getY();
		
		Algebra al = new Algebra();
		
		double up = al.mult(x1, x2) + al.mult(y1, y2);
		double down = al.norm2(x1.assign(y1, Functions.min)) * al.norm2(x2.assign(y2, Functions.min));
		theta = up / down;
		
		return Math.acos(theta);
	}
}
