package bioinfo.proteins;

public enum SecStructEight {

	G ('G',SecStructThree.H),
	H ('H',SecStructThree.H),
	I ('I',SecStructThree.H),
	E ('E',SecStructThree.E),
	B ('B',SecStructThree.E),
	T ('T',SecStructThree.C),
	S ('S',SecStructThree.C),
	C ('C',SecStructThree.C);
	
	private final char charRepres;
	private final SecStructThree threeClassAnalogon;
	private SecStructEight(char charRespres, SecStructThree sst){
		this.charRepres = charRespres;
		this.threeClassAnalogon = sst;
	}
	
	public SecStructThree getThreeClassAnalogon(){
		return threeClassAnalogon;
	}
	
	public char getCharRespres(){
		return charRepres;
	}
	
	public static SecStructEight getSSFromChar(char ss){
		for(SecStructEight ssE: SecStructEight.values()){
			if(ssE.getCharRespres() == ss){
				return ssE;
			}
		}
		return SecStructEight.C;
	}
	
	//addded for readout from MySQL-DB (ask Paul)
	public static SecStructEight getSSFromChar(String ssString){
		char ss;
		if(ssString != null){
			ss = ssString.charAt(0);
		}else{
			return SecStructEight.C;
		}
		
		for(SecStructEight ssE: SecStructEight.values()){
			if(ssE.getCharRespres() == ss){
				return ssE;
			}
		}
		return SecStructEight.C;
	}
}
