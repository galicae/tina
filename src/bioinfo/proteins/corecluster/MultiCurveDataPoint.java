package bioinfo.proteins.corecluster;


public class MultiCurveDataPoint {

	private String id;
	private double voronoiArea;
	private double theta;
	private String type;
	public int minLength;
	
	
	public MultiCurveDataPoint(String id, double voronoiArea, double theta,
			String type, int minLength) {
		super();
		this.id = id;
		this.voronoiArea = voronoiArea;
		this.theta = theta * 180.0 / Math.PI;
		this.type = type;
		this.minLength = minLength;
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
	public int getMinLength() {
		return minLength;
	}
	public void setMinLength(int minLength) {
		this.minLength = minLength;
	}

	

}
