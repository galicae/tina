package test;

//import highscorealignments.HighScoreAlign;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.InputStreamReader;
import java.util.ArrayList;
//import java.util.Map.Entry;

import bioinfo.Sequence;
import bioinfo.alignment.SequenceAlignment;
//import bioinfo.alignment.SequenceAlignmentFileReader;
import bioinfo.alignment.gotoh.Gotoh;
import bioinfo.proteins.PDBEntry;
import bioinfo.proteins.PDBFileReader;
import bioinfo.superpos.TMMain;
import bioinfo.superpos.Transformation;

public class aufg1Test {
	public static void main(String[] args) throws Exception {
		PDBFileReader reader = new PDBFileReader();
		TMMain main = new TMMain();
		FileWriter tmStream = new FileWriter("tmScores");
		FileWriter gdtStream = new FileWriter("gdtScores");
		BufferedWriter tm = new BufferedWriter(tmStream);
		BufferedWriter gdt = new BufferedWriter(gdtStream);
		// SequenceAlignmentFileReader ali = new SequenceAlignmentFileReader();
		PDBEntry pdb1;
		PDBEntry pdb2;
		Transformation tr = new Transformation(null, null, null, null, 0, 0, 0);
		ArrayList<Double> scores = new ArrayList<Double>();

		// // inpairs, outpairs, seqlib
		// HighScoreAlign test = new HighScoreAlign(args[0], args[1], args[2]);
		//
		// // matrix, Gopen, gextend, mode
		// System.out.println("calculating alignments...");
		// test.calculation("dayhoff.mat", -12, -1, "freeshift");
		//
		// System.out.println(test.getHighScoreAligns().entrySet().size());

		// System.out.println("Reading....");
		// ArrayList<SequenceAlignment> alignments = new
		// ArrayList<SequenceAlignment>();
		// alignments =
		// (ArrayList<SequenceAlignment>)ali.readAlignments(args[0]);
		// System.out.println("Writing....");

		BufferedReader br = null;
		String line;
		int count = 0;
		String[] tmp = new String[3];
		Sequence seq1;
		Sequence seq2;
		String seq1tmp = "";
		String seq2tmp = "";
		String ali1tmp;
		String ali2tmp;
		try {
			br = new BufferedReader(new InputStreamReader(new FileInputStream(args[0])));
			while ((line = br.readLine()) != null) {
				tmp[count % 3] = line;
				count++;
				if (count % 3 == 0) {
					ali1tmp = tmp[1].split(":")[1].trim();
					ali2tmp = tmp[2].split(":")[1].trim();
					for (int i = 0; i != ali1tmp.length(); i++) {
						if (ali1tmp.charAt(i) != '-') {
							seq1tmp += ali1tmp.charAt(i);
						}
						if (ali2tmp.charAt(i) != '-') {
							seq2tmp += ali2tmp.charAt(i);
						}
					}
					SequenceAlignment cur = new SequenceAlignment(new Sequence(
							tmp[1].split(":")[0].trim(), seq1tmp),
							new Sequence(tmp[2].split(":")[0].trim(), seq2tmp),
							tmp[1].split(":")[1].trim(),
							tmp[2].split(":")[1].trim(),
							(int) (Gotoh.FACTOR * Double.parseDouble(tmp[0]
									.split("\\s+")[2].trim())));
					pdb1 = reader.readPDBFromFile("C:/Users/nikos/Desktop/STRUCTURES/"
							+ cur.getComponent(0).getId() + ".pdb");
					pdb2 = reader.readPDBFromFile("C:/Users/nikos/Desktop/STRUCTURES/"
							+ cur.getComponent(1).getId() + ".pdb");
					tr = main.calculateTransformation(cur, pdb1, pdb2);
					System.out.println(cur.getComponent(0).getId());
					tm.write(tr.getTmscore() + "\n");
					gdt.write(tr.getGdt() + "\n");
				}
			}
		} catch (Exception e) {
			e.printStackTrace();
		}

		
		tm.close();
		gdt.close();
	}
}
