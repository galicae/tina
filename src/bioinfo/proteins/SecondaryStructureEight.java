package bioinfo.proteins;

public enum SecondaryStructureEight {

	G ('G',SecondaryStructureThree.H),
	H ('H',SecondaryStructureThree.H),
	I ('I',SecondaryStructureThree.H),
	E ('E',SecondaryStructureThree.E),
	B ('B',SecondaryStructureThree.E),
	T ('T',SecondaryStructureThree.C),
	S ('S',SecondaryStructureThree.C),
	C ('C',SecondaryStructureThree.C);
	
	private final char charRepres;
	private final SecondaryStructureThree threeClassAnalogon;
	private SecondaryStructureEight(char charRespres, SecondaryStructureThree sst){
		this.charRepres = charRespres;
		this.threeClassAnalogon = sst;
	}
	
	public SecondaryStructureThree getThreeClassAnalogon(){
		return threeClassAnalogon;
	}
	
	public char getCharRespres(){
		return charRepres;
	}
	
	public static SecondaryStructureEight getSSFromChar(char ss){
		for(SecondaryStructureEight ssE: SecondaryStructureEight.values()){
			if(ssE.getCharRespres() == ss){
				return ssE;
			}
		}
		return SecondaryStructureEight.C;
	}
}
