package test;

import bioinfo.alignment.thr3ad.InitClass;

public class HydropathyTest {
	public static void main(String[] args) {
		InitClass bla = new InitClass();
		
		double[] test = bla.calcHydropathyScores();
		double[][] matrix = bla.calcHydropathyMatrix(test);
		bla.printMatrix(matrix, test);
	}
}
