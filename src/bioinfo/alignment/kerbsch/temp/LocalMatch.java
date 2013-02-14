package bioinfo.alignment.kerbsch.temp;

public class LocalMatch{
	private int score;
	private int[] coord = new int[2];
	
	public LocalMatch(int score,int[] coord){
		this.score = score;
		this.coord = coord;
	}
	
	public int[] getCoords(){
		return this.coord;
	}
	
	public int getScore(){
		return this.score;
	}
	
	public void setScore(int score){
		this.score = score;
	}
}
