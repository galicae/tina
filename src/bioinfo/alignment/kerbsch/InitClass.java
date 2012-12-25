package bioinfo.alignment.kerbsch;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.lang.reflect.Field;

public class InitClass {

	/**
	 * this function reads the scores from HydropathyScores
	 * 
	 * @return the scores from the class (that contains alphabetical values as
	 *         keys and numerics as values; they must all be static fields).The
	 *         matrix is sparse, and can be used for every sequence from the
	 *         latin alphabet
	 */
	public double[] calcHydropathyScores() {
		Field[] vars = HydropathyScores.class.getFields();
		double[] hydropathyScore = new double[26];
		for (Field f : vars) {
			char cur = f.getName().charAt(0);
			try {
				hydropathyScore[cur - 65] = f.getDouble(f);
			} catch (IllegalArgumentException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			} catch (IllegalAccessException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}
		return hydropathyScore;
	}

	/**
	 * this function reads the scores from HydropathyScores
	 * 
	 * @return the scores from the class (that contains alphabetical values as
	 *         keys and numerics as values; they must all be static fields).The
	 *         matrix is sparse, and can be used for every sequence from the
	 *         latin alphabet
	 */
	public double[] calcPolarityScores() {
		Field[] vars = PolarityScores.class.getFields();
		double[] polarityScore = new double[26];
		for (Field f : vars) {
			char cur = f.getName().charAt(0);
			try {
				polarityScore[cur - 65] = f.getDouble(f);
			} catch (IllegalArgumentException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			} catch (IllegalAccessException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}
		return polarityScore;
	}
	
	/**
	 * calculates a matrix that gotoh can accept as a matrix, based on the
	 * hydropathy scores
	 * 
	 * @param hScore
	 *            the hydropathy scores as a sparse double matrix (length 26)
	 */
	public double[][] calcGotohInputMatrix(double[] hScore, int mult, double add) {
		double[][] matrix = new double[hScore.length][hScore.length];
		double temp = 0;
		for (int i = 0; i < hScore.length; i++) {
			if (hScore[i] == 0)
				continue;
			for (int j = 0; j < hScore.length; j++) {
				if (hScore[j] == 0)
					continue;
				temp = hScore[i] - hScore[j];
				temp = (1 - Math.abs(temp / mult)) * mult * 1.0 - add;
				temp = ((int) (temp * 1000)) / 1000.0;
				matrix[i][j] = temp;
			}
		}
		return matrix;
	}

	public void printMatrix(double[][] matrix, double[] hScore, String name) {
		try {
			BufferedWriter writer = new BufferedWriter(new FileWriter(
					"hydro.mat"));
			writer.write("NAME " + name + "\n");
			writer.write("SCORE similarity\n");
			writer.write("NUMROW dd" + "\n");
			writer.write("NUMCOL dd" + "\n");
			writer.write("ROWINDEX ");
			for (int i = 0; i < hScore.length; i++) {
				if (matrix[0][i] != 0) {
					writer.write((char) (65 + i));
				}
			}

			writer.write("\nCOLINDEX ");
			for (int i = 0; i < hScore.length; i++) {
				if (matrix[i][0] != 0) {
					writer.write((char) (65 + i));
				}
			}
			writer.write("\n");

			for (int i = 0; i < hScore.length; i++) {
				if (hScore[i] == 0)
					continue;
				writer.write("MATRIX \t");
				for (int j = 0; j < hScore.length; j++) {
					if (hScore[j] == 0)
						continue;
					writer.write(matrix[j][i] + "\t");
				}
				writer.write("\n");
			}
			writer.close();
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}
