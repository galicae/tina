package bioinfo.alignment.kerbsch.temp;

import java.util.ArrayList;
import java.util.List;

public class LocalCore {
	private int score;
	private List<int[]> coords = new ArrayList<int[]>();
	
	public LocalCore(int score, List<int[]> coords){
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
}
