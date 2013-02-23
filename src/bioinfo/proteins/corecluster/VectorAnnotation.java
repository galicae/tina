package bioinfo.proteins.corecluster;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.List;

public class VectorAnnotation implements Serializable{

	/**
	 * 
	 */
	private static final long serialVersionUID = -9038314963462330443L;
	private String id;
	private List<Curve> curves;

	public VectorAnnotation(String id) {
		this.curves = new ArrayList<Curve>();
		this.id = id;
	}

	public void appendVector(Curve curve) {
		this.curves.add(curve);
	}

	public void appendVectors(List<Curve> curves) {
		this.curves.addAll(curves);
	}

	public List<Curve> getAllVectors() {
		return curves;
	}

	public String getId() {
		return id;
	}
}
