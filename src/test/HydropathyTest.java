package test;

import bioinfo.alignment.kerbsch.InitClass;

public class HydropathyTest {
	public static void main(String[] args) throws ClassNotFoundException {
		InitClass bla = new InitClass();
		double[] test = bla.calcPolarityScores();
		double[][] matrix = bla.calcGotohInputMatrix(test);
//		System.out.println(matrix.length);
		bla.printMatrix(matrix, test, "polarity");
	}
}
