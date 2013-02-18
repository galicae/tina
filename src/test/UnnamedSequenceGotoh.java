package test;

import bioinfo.Sequence;
import bioinfo.alignment.SequenceAlignment;
import bioinfo.alignment.gotoh.FreeshiftSequenceGotoh;
import bioinfo.alignment.gotoh.Gotoh;
import bioinfo.alignment.gotoh.LocalSequenceGotoh;
import bioinfo.alignment.matrices.QuasarMatrix;

public class UnnamedSequenceGotoh {

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		String id1 = args[0];
		String id2 = args[2];
		String seq1 = args[1];
		String seq2 = args[3];
		Sequence sequence1 = new Sequence(id1, seq1);
		Sequence sequence2 = new Sequence(id2, seq2);
		
		double[][] matrix = QuasarMatrix.DAYHOFF_MATRIX;
		FreeshiftSequenceGotoh gotoh = new FreeshiftSequenceGotoh(-12, -1, matrix);
		SequenceAlignment ali = gotoh.align(sequence1, sequence2);
		System.out.print(String.format("%.3f",ali.getScore()));
	}

}
