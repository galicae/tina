package test;

import java.util.HashMap;

import bioinfo.Sequence;
import bioinfo.alignment.SequenceAlignment;
import bioinfo.alignment.gotoh.FreeshiftSequenceGotoh;
import bioinfo.alignment.kerbsch.temp.SeqLibrary;
import bioinfo.alignment.matrices.QuasarMatrix;
import bioinfo.proteins.PDBEntry;
import bioinfo.proteins.PDBFileReader;
import bioinfo.superpos.TMMain;
import bioinfo.superpos.Transformation;

public class TMPipeline {
	public static void main(String[] args) throws Exception {

		// 1wjfA00 1zr4B03

		String seq1 = "1lddB00";
		String seq2 = "1wq2B00";
		// String seq2 = "1lddB00";

		HashMap<String, char[]> sequences = SeqLibrary.read("domains.seqlib");
		double[][] matr = QuasarMatrix.parseMatrix("dayhoff.mat");

		PDBFileReader reader = new PDBFileReader();
		PDBEntry pdb1 = reader.readPDBFromFile("./" + seq1 + ".pdb");
		PDBEntry pdb2 = reader.readPDBFromFile("./" + seq2 + ".pdb");

		// System.out.println(pdb1.getAtomSectionAsString());
		FreeshiftSequenceGotoh gotoh = new FreeshiftSequenceGotoh(-12, -1, matr);
		Sequence sequence1 = new Sequence(seq1, sequences.get(seq1));
		Sequence sequence2 = new Sequence(seq2, sequences.get(seq2));

		SequenceAlignment aliSeq = gotoh.align(sequence1, sequence2);

		System.out.println(aliSeq.toStringVerbose());
		double cur = 0;
		double prev = 0;
		boolean tmAlright = true;
		TMMain main = new TMMain();
		Transformation tr1 = main.calculateTransformation(aliSeq, pdb1,
				pdb2);
		prev = tr1.getTmscore();
		for (int i = 0; i < 10000; i++) {
			main = new TMMain();
			tr1 = main.calculateTransformation(aliSeq, pdb1,
					pdb2);
			cur = tr1.getTmscore();
			if(cur != prev) {
				tmAlright = false;
			}
			prev = cur;
		}
		
		if(tmAlright)
			System.out.println("tm score is alright");
		else
			System.out.println("BEEP BEEP");
	}
}
