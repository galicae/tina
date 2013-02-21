package test;


import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.InputStreamReader;
import java.util.HashMap;

import bioinfo.Sequence;
import bioinfo.alignment.SequenceAlignment;
//import bioinfo.alignment.SequenceAlignmentFileReader;
import bioinfo.alignment.gotoh.FreeshiftSequenceGotoh;
import bioinfo.alignment.kerbsch.temp.SeqLibrary;
import bioinfo.alignment.matrices.QuasarMatrix;
import bioinfo.proteins.PDBEntry;
import bioinfo.proteins.PDBFileReader;
import bioinfo.superpos.TMMain;
import bioinfo.superpos.Transformation;

public class aufg1Test2 {
	public static void main(String[] args) throws Exception {
		FileWriter tmStream = new FileWriter("tmScores");
		FileWriter gdtStream = new FileWriter("gdtScores");
		BufferedWriter tm = new BufferedWriter(tmStream);
		BufferedWriter gdt = new BufferedWriter(gdtStream);
		// SequenceAlignmentFileReader ali = new SequenceAlignmentFileReader();
		
		BufferedReader br = new BufferedReader(new InputStreamReader(new FileInputStream(args[1])));
		String line = null;
		
		double[][] matrix = QuasarMatrix.parseMatrix(args[2]);
		
		HashMap<String, char[]> sequences = SeqLibrary.read(args[0]);
		
		while ((line = br.readLine()) != null) {
			if(line.startsWith(">")) {
				String[] ids = line.substring(1).split("\\s+");
				FreeshiftSequenceGotoh gotoh = new FreeshiftSequenceGotoh(-12, -1, matrix);
				Sequence sequence1 = new Sequence(ids[0], sequences.get(ids[0]));
				Sequence sequence2 = new Sequence(ids[1], sequences.get(ids[1]));
				
				SequenceAlignment ali = gotoh.align(sequence1, sequence2);
				
				PDBFileReader reader = new PDBFileReader();
				PDBEntry pdb1 = reader.readPDBFromFile("C:/Users/nikos/Desktop/STRUCTURES/" + ali.getComponent(0).getId() + ".pdb");
				PDBEntry pdb2 = reader.readPDBFromFile("C:/Users/nikos/Desktop/STRUCTURES/" + ali.getComponent(1).getId() + ".pdb");
				
				TMMain main = new TMMain();
				Transformation tr = main.calculateTransformation(ali, pdb1, pdb2);
				System.out.println(ali.getComponent(0).getId() + " " + ali.getComponent(1).getId() + " " + tr.getTmscore());
//				System.out.println(tr.getTmscore());
				tm.write(tr.getTmscore() + "\n");
				gdt.write(tr.getGdt() + "\n");
				if(tr.getTmscore() >= 1) {
					System.err.println(ali.getComponent(0).getId() + " " + ali.getComponent(1).getId() + " " + tr.getTmscore());
					break;
				}
			}
		}
		br.close();
		tm.close();
		gdt.close();
	}
}
