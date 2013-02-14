package bioinfo.energy.potential.voroEval;

import bioinfo.proteins.AminoAcidName;

public class VoroEvalDataPoint {
	
	private AminoAcidName amino;
	private String sourceId;
	private int sourceIndex;
	private double dsspAcc;
	private double selfAcc;
	
	public VoroEvalDataPoint(AminoAcidName amino, String sourceId, int sourceIndex, double dsspAcc, double selfAcc){
		this.amino = amino;
		this.sourceId = sourceId;
		this.sourceIndex = sourceIndex;
		this.dsspAcc = dsspAcc;
		this.selfAcc = selfAcc;
	}

	public AminoAcidName getAmino() {
		return amino;
	}

	public void setAmino(AminoAcidName amino) {
		this.amino = amino;
	}

	public String getSourceId() {
		return sourceId;
	}

	public void setSourceId(String sourceId) {
		this.sourceId = sourceId;
	}

	public int getSourceIndex() {
		return sourceIndex;
	}

	public void setSourceIndex(int sourceIndex) {
		this.sourceIndex = sourceIndex;
	}

	public double getDsspAcc() {
		return dsspAcc;
	}

	public void setDsspAcc(double dsspAcc) {
		this.dsspAcc = dsspAcc;
	}

	public double getSelfAcc() {
		return selfAcc;
	}

	public void setSelfAcc(double selfAcc) {
		this.selfAcc = selfAcc;
	}
	
	

}
