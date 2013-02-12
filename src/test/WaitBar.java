package test;

import java.util.ArrayList;
import java.util.List;

import bioinfo.Sequence;
import bioinfo.alignment.SequenceAlignment;
import bioinfo.alignment.gotoh.GlobalSequenceGotoh;
import bioinfo.alignment.matrices.QuasarMatrix;
import bioinfo.proteins.PDBEntry;
import bioinfo.proteins.PDBFileReader;

public class WaitBar {
	public static void main(String[] args) {
		PDBFileReader reader = new PDBFileReader("./proteins2");
		List<PDBEntry> list = new ArrayList<PDBEntry>();
		ArrayList<SequenceAlignment> alignments = new ArrayList<SequenceAlignment>();
		list = reader.readPdbFolder();
		double[][] matrix = QuasarMatrix.DAYHOFF_MATRIX;

		for (int i = 0; i < list.size(); i++) {
			PDBEntry p = list.get(i);
			String seq = "";
			for (int k = 0; k < p.length(); k++) {
				seq += p.getAminoAcid(i);
			}
			Sequence seq1 = new Sequence(p.getID(), seq);
			for (int j = 0; j < list.size(); j++) {
				if(j == i)
					continue;
				PDBEntry q = list.get(j);
				seq = "";
				for (int k = 0; k < q.length(); k++) {
					seq += q.getAminoAcid(i);
				}
				Sequence seq2 = new Sequence(q.getID(), seq);

				GlobalSequenceGotoh gotoh = new GlobalSequenceGotoh(-9, -1,
						matrix);
				alignments.add(gotoh.align(seq1, seq2));
			}
		}

		// now evaluate alignments for sequence identity
		double match = 0;
		double total = 0;
		for(SequenceAlignment a: alignments) {
			for(int i = 0; i < a.length(); i++) {
				if(a.getRow(0)[i] == a.getRow(1)[i])
					match++;
				total++;
			}
		}
		System.out.println(match / total);
	}
}