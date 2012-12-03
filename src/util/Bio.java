package util;

import java.util.HashMap;

public class Bio {
	
	public static HashMap<String,Character> threeLetterCode = new HashMap<String,Character>();
	public static HashMap<Character,String> oneLetterCode = new HashMap<Character,String>();
	
	static{
		threeLetterCode.put("ALA",'A');
		threeLetterCode.put("ARG",'R');
		threeLetterCode.put("ASN",'N');
		threeLetterCode.put("ASP",'D');
		threeLetterCode.put("CYS",'C');
		threeLetterCode.put("GLU",'E');
		threeLetterCode.put("GLN",'Q');
		threeLetterCode.put("GLY",'G');
		threeLetterCode.put("HIS",'H');
		threeLetterCode.put("ILE",'I');
		threeLetterCode.put("LEU",'L');
		threeLetterCode.put("LYS",'K');
		threeLetterCode.put("MET",'M');
		threeLetterCode.put("PHE",'F');
		threeLetterCode.put("PRO",'P');
		threeLetterCode.put("SER",'S');
		threeLetterCode.put("THR",'T');
		threeLetterCode.put("TRP",'W');
		threeLetterCode.put("TYR",'Y');
		threeLetterCode.put("VAL",'V');
		threeLetterCode.put("UNK",'X');
		
		
		oneLetterCode.put('A',"ALA");
		oneLetterCode.put('R',"ARG");
		oneLetterCode.put('N',"ASN");
		oneLetterCode.put('D',"ASP");
		oneLetterCode.put('C',"CYS");
		oneLetterCode.put('E',"GLU");
		oneLetterCode.put('Q',"GLN");
		oneLetterCode.put('G',"GLY");
		oneLetterCode.put('H',"HIS");
		oneLetterCode.put('I',"ILE");
		oneLetterCode.put('L',"LEU");
		oneLetterCode.put('K',"LYS");
		oneLetterCode.put('M',"MET");
		oneLetterCode.put('F',"PHE");
		oneLetterCode.put('P',"PRO");
		oneLetterCode.put('S',"SER");
		oneLetterCode.put('T',"THR");
		oneLetterCode.put('W',"TRP");
		oneLetterCode.put('Y',"TYR");
		oneLetterCode.put('V',"VAL");
		oneLetterCode.put('X',"UNK");



	}
		

	
	public static Character codeTranslate(String in){
		return threeLetterCode.get(in);
	}
	
	public static String codeTranslate(Character in){
		return oneLetterCode.get(in);
	}

}
