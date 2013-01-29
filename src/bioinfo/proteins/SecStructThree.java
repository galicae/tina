package bioinfo.proteins;

public enum SecStructThree {

	H ('H'),
	E ('E'),
	C ('C');
	
	private final char charRepres;
	private SecStructThree(char charRepres){
		this.charRepres = charRepres;
	}
	
	public char getCharRepres(){
		return this.charRepres;
	}
	
}
