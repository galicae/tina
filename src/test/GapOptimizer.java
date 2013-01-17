package test;


import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.InputStreamReader;
import java.util.HashMap;

import bioinfo.Sequence;
import bioinfo.alignment.SequenceAlignment;
import bioinfo.alignment.gotoh.FreeshiftSequenceGotoh;
import bioinfo.alignment.kerbsch.SeqLibrary;
import bioinfo.alignment.matrices.QuasarMatrix;
import bioinfo.proteins.PDBEntry;
import bioinfo.proteins.PDBFileReader;
import bioinfo.superpos.TMMain;
import bioinfo.superpos.Transformation;

/**
 * this class is meant to calculate the optimal gap conditions for a given
 * matrix
 * 
 * @author papadopoulos
 * 
 */
public class GapOptimizer {

	/**
	 * the main function here takes a seqlib file at 0, the pairfile in 1, the
	 * matrix in question in 2 and the gap parameters at 3 and 4. It is supposed
	 * to be integrated in a bash script.
	 * 
	 * @param args
	 *            0: seqlib file, 1: pairfile, 2: gotoh matrix, 3: gap open, 4: gap extend
	 * @throws Exception
	 */
	public static void main(String[] args) throws Exception {

		HashMap<String, char[]> sequences = SeqLibrary.read(args[0]);
		double[][] hydMat = QuasarMatrix.parseMatrix(args[2]);
		int go = Integer.parseInt(args[3]);
		int ge = Integer.parseInt(args[4]);

		int seq = 0;
		long scoreSum = 0;

		TMMain tmCalculator = new TMMain();
		Transformation tr = new Transformation(null, null, null, null, 0, 0, 0);
		
		BufferedReader br = new BufferedReader(new InputStreamReader(
				new FileInputStream(args[1])));
		String line = null;
		while ((line = br.readLine()) != null) {
			String[] ids = line.split("\\s+");
			seq++;
			FreeshiftSequenceGotoh gotohHyd = new FreeshiftSequenceGotoh(go,
					ge, hydMat);
			Sequence sequence1 = new Sequence(ids[0], sequences.get(ids[0]));
			Sequence sequence2 = new Sequence(ids[1], sequences.get(ids[1]));
			
			SequenceAlignment aliHyd = gotohHyd.align(sequence1, sequence2);
			
			PDBFileReader reader = new PDBFileReader();
			PDBEntry pdb1 = reader.readPDBFromFile("/home/p/papadopoulos/Desktop/STRUCTURES/" + aliHyd.getComponent(0).getID() + ".pdb");
			PDBEntry pdb2 = reader.readPDBFromFile("/home/p/papadopoulos/Desktop/STRUCTURES/" + aliHyd.getComponent(1).getID() + ".pdb");
			
			tr= tmCalculator.calculateTransformation(aliHyd, pdb1, pdb2);
			System.out.println(tr.getTmscore());
		}
		br.close();
		
	}
}
