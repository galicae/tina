package bioinfo.alignment.kerbsch.temp;

import java.util.ArrayList;
import java.util.List;

public class Locals {
	private int score;
	private double evalue;
	private List<int[]> coords = new ArrayList<int[]>();
	
	public Locals(int score, List<int[]> coords){
		this.score = score;
		this.coords = coords;
	}
	
	public List<int[]> getCoords(){
		return this.coords;
	}
	
	public int getScore(){
		return this.score;
	}
	
	public void setScore(int score){
		this.score = score;
	}
	
	public void setEvalue(double evalue){
		this.evalue = evalue;
	}
	
	public double getEvalue(){
		return this.evalue;
	}
}
