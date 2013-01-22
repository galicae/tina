package test;

import bioinfo.alignment.kerbsch.temp.InitClass;

public class HydropathyTest {
	public static void main(String[] args) throws ClassNotFoundException {
		InitClass bla = new InitClass();
		double[] test = bla.calcSecStructScores();
		double[][] matrix = bla.calcGotohInputMatrix(test);
//		System.out.println(matrix.length);
		bla.printMatrix(matrix, test, "secstruct");
	}
}
