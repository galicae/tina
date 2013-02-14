package bioinfo.proteins.fragm3nt.run;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.util.LinkedList;

import bioinfo.Sequence;
import bioinfo.alignment.SequenceAlignment;
import bioinfo.alignment.gotoh.FreeshiftSequenceGotoh;
import bioinfo.alignment.matrices.QuasarMatrix;

public class AssessData {

	public static void main(String[] args) {
		try {
			// read pairs
			LinkedList<String[]> pairs = new LinkedList<String[]>();
			readPairs(pairs);
			System.out.println("read pairs");
			// read sequences
			LinkedList<String> seqs = new LinkedList<String>();
			readSequences(seqs);
			System.out.println("read sequences");
			BufferedWriter list09 = new BufferedWriter(new FileWriter("bucket09"));
			BufferedWriter list07 = new BufferedWriter(new FileWriter("bucket07"));
			BufferedWriter list05 = new BufferedWriter(new FileWriter("bucket05"));
			double[] scores = new double[pairs.size()];

			for (int i = 0; i < pairs.size(); i++) {
				System.out.println(i + " of " + pairs.size());
				scores[i] = align(pairs.get(i)[0], pairs.get(i)[1], seqs);
				if (scores[i] >= 0.5) {
					list05.write(pairs.get(i)[0] + " " + pairs.get(i)[1] + " " + scores[i] + "\n");
					if (scores[i] >= 0.7) {
						list07.write(pairs.get(i)[0] + " " + pairs.get(i)[1] + " " + scores[i] + "\n");
						if (scores[i] >= 0.9) {
							list09.write(pairs.get(i)[0] + " " + pairs.get(i)[1] + " " + scores[i] + "\n");
						}
					}
				}
			}
			list05.close();
			list07.close();
			list09.close();
		} catch (Exception e) {
			e.printStackTrace();
		}
	}

	public static void readPairs(LinkedList<String[]> pairs) {
		try {
			BufferedReader r = new BufferedReader(new FileReader(
					"cathscop.inpairs"));
			String line = "";
			String[] arr = new String[2];
			while ((line = r.readLine()) != null) {
				arr = line.split(" ");
				String[] p = { arr[0], arr[1] };
				pairs.add(p);
			}
			r.close();
		} catch (Exception e) {
			e.printStackTrace();
		}
	}

	public static void readSequences(LinkedList<String> seqs) {
		try {
			BufferedReader r = new BufferedReader(new FileReader(
					"domains.seqlib"));
			String line;
			while ((line = r.readLine()) != null) {
				seqs.add(line);
			}
			r.close();
		} catch (Exception e) {
			e.printStackTrace();
		}
	}

	/**
	 * this function calculates the alignment and returns the sequence identity
	 * that results from it
	 * 
	 * @param id1
	 *            ID of seq1
	 * @param id2
	 *            ID of seq2
	 * @param seqs
	 *            a list containing all sequences in cathscop.inpairs
	 * @return
	 */
	public static double align(String id1, String id2, LinkedList<String> seqs) {
		String seq1 = "";
		String seq2 = "";

		for (String s : seqs) {
			if (s.startsWith(id1))
				seq1 = s.split(":")[1];
			if (s.startsWith(id2))
				seq2 = s.split(":")[1];
		}

		Sequence sequ1 = new Sequence(id1, seq1.toCharArray());
		Sequence sequ2 = new Sequence(id2, seq2.toCharArray());
		FreeshiftSequenceGotoh fg = new FreeshiftSequenceGotoh(-9, -1,
				QuasarMatrix.DAYHOFF_MATRIX);

		SequenceAlignment ali = fg.align(sequ1, sequ2);
		char[][] arr = new char[2][];
		arr[0] = ali.getRow(0);
		arr[1] = ali.getRow(1);
		int match = 0;
		for (int i = 0; i < arr[0].length; i++) {
			if (arr[0][i] == arr[1][i])
				match++;
		}

		double result = (seq1.length() + seq2.length()) / 2;
		result = match / result;
		return result;
	}
}
