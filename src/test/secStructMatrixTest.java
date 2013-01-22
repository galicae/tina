package test;

import bioinfo.alignment.kerbsch.temp.SecStructScores;

public class secStructMatrixTest {

	public static void main(String[] args) {
		System.out.print("\t");
		for (int i = 0; i < SecStructScores.matrix.length; i++) {
			if((char)(i+65)=='C' || (char)(i+65)=='E' || (char)(i+65)=='H'){
				System.out.print((char)(i+65) + "\t");
			}
		}
		
		System.out.println();
				
		for (int i = 0; i < SecStructScores.matrix[0].length; i++) {
			if((char)(i+65)=='C' || (char)(i+65)=='E' || (char)(i+65)=='H'){
				System.out.print((char)(i+65) + "\t");
			}
			for (int j = 0; j < SecStructScores.matrix.length; j++) {
				if(SecStructScores.matrix[i][j] > -999.0){
					System.out.print(SecStructScores.matrix[i][j] + "\t");
				}
			}
			if((char)(i+65)=='C' || (char)(i+65)=='E' || (char)(i+65)=='H'){
				System.out.println();
			}
		}
	}

}
