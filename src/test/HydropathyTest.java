package test;

import bioinfo.alignment.kerbsch.InitClass;

public class HydropathyTest {
	public static void main(String[] args) throws ClassNotFoundException {
		InitClass bla = new InitClass();
		double[] test = bla.calcPolarityScores();
		double min = 0;
		double max = 0;
		for (int i = 0; i < test.length; i++) {
			if (min > test[i])
				min = test[i];
			if (max < test[i])
				max = test[i];
		}
		int mult = (int) (max - min + 1);
		double add = (max - min) / 2.0;
		double[][] matrix = bla.calcGotohInputMatrix(test, mult, add);
		bla.printMatrix(matrix, test, "polarity");

	}
}
