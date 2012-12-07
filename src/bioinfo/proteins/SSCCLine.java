package bioinfo.proteins;

public class SSCCLine {
	
	private final AminoAcidName name;
	private char secStruct;
	private int localContacts;
	private int globalContacts;

	public SSCCLine(AminoAcidName name, char secStruct, int locCont, int globCont){
		this.name = name;
		this.secStruct = secStruct;
		this.localContacts = locCont;
		this.globalContacts = globCont;
	}
	
	public AminoAcidName getAmino(){
		return this.name;
	}
	
	public char getSecStruct(){
		return this.secStruct;
	}
	
	public int getGlobCont(){
		return this.globalContacts;
	}
	
	public int getLocCont(){
		return this.localContacts;
	}
}
