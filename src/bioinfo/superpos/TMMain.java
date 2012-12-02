package bioinfo.superpos;
import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;

import EDU.oswego.cs.dl.util.concurrent.Rendezvous.Rotator;
import bioinfo.proteins.AminoAcid;
import bioinfo.proteins.PDBEntry;
import cern.colt.matrix.DoubleFactory1D;
import cern.colt.matrix.DoubleFactory2D;
import cern.colt.matrix.DoubleMatrix1D;
import cern.colt.matrix.DoubleMatrix2D;

public class TMMain {

	public void calculateTMScore(String alignment, String pFile, String qFile)
			throws Exception {
		TMCollective main = new TMCollective();
		PDBEntry[] pdbs = main.createTMInput(alignment, pFile, qFile);

		writeToFile("TM" + pFile, pdbs[0]);
		writeToFile("TM" + qFile, pdbs[1]);

		// actual calculation of TM score and corresponding rotation matrix
		double[][] tmResult = TMScore.doStuff("TM" + pFile + " TM" + qFile);

		// remember that rmResult is [5][4], and that [i][0] is empty
		// also [4][0] is the TM score and [4][1] the GDT
		// read the matrices and scores:
		DoubleMatrix2D R = readRotationMatrix(tmResult);
		DoubleMatrix1D T = readTMatrix(tmResult);
		
//		Rotator rot = new Rotator(R, T, qFile);
		Transformation tr = new Transformation(null, null, T, R, tmResult[4][2], tmResult[4][1], tmResult[4][0]);
		tr.transform(pdbs[1]);
	}

	/**
	 * read the first row of the multi-result double[][] I get from TMScore. It
	 * contains in [0][i] the T vector. The [0][0] position is empty because
	 * Java vectors start at 1 in the Zhang lab.
	 * 
	 * @param tmResult
	 *            a 5x4 double array containing the R, T matrices as well as the
	 *            TM and GDT scores
	 * @return the T vector
	 */
	private DoubleMatrix1D readTMatrix(double[][] tmResult) {
		double[] doubleT = new double[3];
		for (int i = 1; i < 4; i++) {
			doubleT[i - 1] = tmResult[0][i];
		}
		DoubleFactory1D factory = DoubleFactory1D.dense;
		DoubleMatrix1D T = factory.make(doubleT);
		return T;
	}

	/**
	 * Yet another function aiming to decongest the main method. This one prints
	 * a TMScore input file.
	 * 
	 * @param name
	 *            the name of the TMScore input file. It gets a TM prefix to
	 *            distinguish it from the "normal" PDB file.
	 * @param input
	 *            the amino acid array. Contains the coordinates of interest
	 * @throws IOException
	 */
	public static void writeToFile(String name, PDBEntry input)
			throws IOException {
		FileWriter outStream = new FileWriter(name);
		BufferedWriter out = new BufferedWriter(outStream);

		// write to file
		for (int i = 0; i < input.length(); i++) {
			out.write(input.getAminoAcid(i).toString() + "\n");
		}
		out.close();
	}

	/**
	 * read rows 2,3 and 4 from the multi-result double[][] I get from TMScore.
	 * They contain the rotation matrix. The [i][0] positions are empty because
	 * Java vectors start at 1 in the Zhang lab.
	 * 
	 * @param tmResult
	 *            a 5x4 double array containing the R, T matrices as well as the
	 *            TM and GDT scores
	 * @return the rotation matrix
	 */
	public static DoubleMatrix2D readRotationMatrix(double[][] tmResult) {
		double[][] doubleR = new double[3][3];
		for (int i = 1; i < 4; i++) {
			for (int j = 1; j < 4; j++) {
				doubleR[i - 1][j - 1] = tmResult[i][j];
			}
		}
		DoubleFactory2D factory = DoubleFactory2D.dense;
		DoubleMatrix2D R = factory.make(doubleR);
		return R;
	}
}
