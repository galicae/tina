package test;


import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.InputStreamReader;
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

public class HydroTest2 {
	public static void main(String[] args) throws Exception {
		// SequenceAlignmentFileReader ali = new SequenceAlignmentFileReader();
		BufferedWriter out = new BufferedWriter(new FileWriter("seqVShyq"));
		TMMain tmCalculator = new TMMain();
		Transformation trSeq = new Transformation(null, null, null, null, 0, 0, 0);
		Transformation trHyd = new Transformation(null, null, null, null, 0, 0, 0);
		
		BufferedReader br = new BufferedReader(new InputStreamReader(new FileInputStream(args[1])));
		String line = null;
		
		double[][] matrix = QuasarMatrix.parseMatrix(args[2]);
		double[][] hydMat = QuasarMatrix.parseMatrix(args[3]);
		int seq = 0;
		int hyd = 0;
		HashMap<String, char[]> sequences = SeqLibrary.read(args[0]);
		
		while ((line = br.readLine()) != null) {
//			if(line.startsWith(">")) {
				String[] ids = line.split("\\s+");
				FreeshiftSequenceGotoh gotohSeq = new FreeshiftSequenceGotoh(-12, -1, matrix);
				FreeshiftSequenceGotoh gotohHyd = new FreeshiftSequenceGotoh(-6, -4, hydMat);
				Sequence sequence1 = new Sequence(ids[0], sequences.get(ids[0]));
				Sequence sequence2 = new Sequence(ids[1], sequences.get(ids[1]));
				
				SequenceAlignment aliSeq = gotohSeq.align(sequence1, sequence2);
				SequenceAlignment aliHyd = gotohHyd.align(sequence1, sequence2);
				
				PDBFileReader reader = new PDBFileReader();
				PDBEntry pdb1 = reader.readPDBFromFile("/home/p/papadopoulos/Desktop/STRUCTURES/" + aliSeq.getComponent(0).getId() + ".pdb");
				PDBEntry pdb2 = reader.readPDBFromFile("/home/p/papadopoulos/Desktop/STRUCTURES/" + aliSeq.getComponent(1).getId() + ".pdb");
				
				trSeq = tmCalculator.calculateTransformation(aliSeq, pdb1, pdb2);
				trHyd = tmCalculator.calculateTransformation(aliHyd, pdb1, pdb2);
//				System.out.print(aliSeq.getComponent(0).getID() + " " + aliSeq.getComponent(1).getID() + " ");
//				System.out.print(trSeq.getTmscore() + "\t" + trHyd.getTmscore() + "\n");
				System.out.println(aliSeq.toStringVerbose());
				System.out.println(aliHyd.toStringVerbose());
				if(trSeq.getTmscore() > trHyd.getTmscore())
					seq++;
				else
					hyd++;
//			}
		}
		br.close();
		out.flush();
		out.close();
		System.err.println(seq + " " + hyd);
	}
}
