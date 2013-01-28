package bioinfo.proteins;

public enum SecondaryStructureThree {

	H ('H'),
	E ('E'),
	C ('C');
	
	private final char charRepres;
	private SecondaryStructureThree(char charRepres){
		this.charRepres = charRepres;
	}
	
	public char getCharRepres(){
		return this.charRepres;
	}
	
}
