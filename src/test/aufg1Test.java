package test;

import highscorealignments.HighScoreAlign;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.util.ArrayList;
import java.util.Map.Entry;

import bioinfo.alignment.SequenceAlignment;
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
		PDBEntry pdb1;
		PDBEntry pdb2;
		Transformation tr = new Transformation(null, null, null, null, 0, 0, 0);
		ArrayList<Double> scores = new ArrayList<Double>();

		// inpairs, outpairs, seqlib
		HighScoreAlign test = new HighScoreAlign(args[0], args[1], args[2]);
		
		// matrix, Gopen, gextend, mode
		System.out.println("calculating alignments...");
		test.calculation("dayhoff.mat", -12, -1, "freeshift");
		
		// the hashmap: String id is the key. Iteration is there!
		System.out.println("calculating and writing...");
		for (Entry<String, SequenceAlignment> a : test.getHighScoreAligns()
				.entrySet()) {
			pdb1 = reader.readPDBFromFile("C:/Users/nikos/Desktop/STRUCTURES/"
					+ a.getValue().getComponent(0).getID() + ".pdb");
			pdb2 = reader.readPDBFromFile("C:/Users/nikos/Desktop/STRUCTURES/"
					+ a.getValue().getComponent(1).getID() + ".pdb");
			SequenceAlignment al = a.getValue();
			tr = main.calculateTransformation(al, pdb1, pdb2);
			tm.write(tr.getTmscore() + "\n");
			gdt.write(tr.getGdt() + "\n");
		}
		tm.close();
		gdt.close();
	}
}
