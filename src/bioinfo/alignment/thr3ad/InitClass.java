package bioinfo.alignment.thr3ad;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.lang.reflect.Field;

public class InitClass {

	/**
	 * this function reads the scores from HydropathyScores
	 * 
	 * @return the scores from HydropathyScores as a double[]; the matrix is
	 *         sparse, and can be used for every sequence from the latin
	 *         alphabet
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
	 * calculates a matrix that gotoh can accept as a matrix, based on the
	 * hydropathy scores
	 * 
	 * @param hScore
	 *            the hydropathy scores as a sparse double matrix (length 26)
	 */
	public double[][] calcHydropathyMatrix(double[] hScore) {
		double[][] matrix = new double[hScore.length][hScore.length];
		double temp = 0;
		for (int i = 0; i < hScore.length; i++) {
			if (hScore[i] == 0)
				continue;
			for (int j = 0; j < hScore.length; j++) {
				if (hScore[j] == 0)
					continue;
				temp = hScore[i] - hScore[j];
				temp = (1 - Math.abs(temp / 9.0)) * 9 - 4.5;
				temp = ((int) (temp * 1000)) / 1000.0;
				matrix[i][j] = temp;
			}
		}
		return matrix;
	}

	public void printMatrix(double[][] matrix, double[] hScore) {
		try {
			BufferedWriter writer = new BufferedWriter(new FileWriter(
					"hydro.mat"));
			writer.write("NAME hydropathy\n");
			writer.write("SCORE similarity\n");
			writer.write("NUMROW 20" + "\n");
			writer.write("NUMCOL 20" + "\n");
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
