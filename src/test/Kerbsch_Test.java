package test;

import bioinfo.alignment.kerbsch.InitClass;
import bioinfo.alignment.kerbsch.Kerbsch;
import bioinfo.alignment.matrices.QuasarMatrix;

public class Kerbsch_Test {

	public static void main(String[] args) {
		InitClass matrices = new InitClass();
		double[][] polMatrix = matrices.calcGotohInputMatrix(matrices.calcPolarityScores());
//		double[][] hbMatrix  = matrices.calcGotohInputMatrix(matrices.calcHydropathyScores());
		double[][] substMatrix = QuasarMatrix.parseMatrix(args[0]);
		for (int i = 0;  i< polMatrix.length; i++) {
			for (int j = 0; j < polMatrix[0].length; j++) {
				System.out.print(polMatrix[i][j] + "\t");
			}
			System.out.println();
		}
		//Kerbsch test = new Kerbsch(-12.0, -1.0, substMatrix, hbMatrix, polMatrix, secStructMatrix);

	}

}
