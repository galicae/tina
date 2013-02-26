package bioinfo.proteins.corecluster;


public class MultiCurveDataPoint {

	private String id;
	private double voronoiArea;
	private double theta;
	private String type;

	public MultiCurveDataPoint(String id, double voronoiArea, double theta,
			String type) {
		this.id = id;
		this.voronoiArea = voronoiArea;
		this.theta = theta / Math.PI * 180;
		this.type = type;
	}

	public String getId() {
		return id;
	}

	public void setId(String id) {
		this.id = id;
	}

	public double getVoronoiArea() {
		return voronoiArea;
	}

	public void setVoronoiArea(double voronoiArea) {
		this.voronoiArea = voronoiArea;
	}

	public double getTheta() {
		return theta;
	}

	public void setTheta(double theta) {
		this.theta = theta;
	}

	public String getType() {
		return type;
	}

	public void setType(String type) {
		this.type = type;
	}

}
