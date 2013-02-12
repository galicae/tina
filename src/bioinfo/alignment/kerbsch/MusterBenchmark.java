package bioinfo.alignment.kerbsch;

import highscorealignments.CathScopEntry;
import highscorealignments.CathScopHash;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;

import bioinfo.alignment.Aligner;
import bioinfo.alignment.kerbsch.temp.InitClass;
import bioinfo.alignment.kerbsch.temp.PairReader;
import bioinfo.alignment.kerbsch.temp.SecStructScores;
import bioinfo.alignment.kerbsch.temp.SeqLibrary;
import bioinfo.alignment.matrices.QuasarMatrix;


public class MusterBenchmark {
	public static void main(String[] args) throws IOException {
		InitClass scores = new InitClass();
		Aligner gotoh = null;
		
		HashMap<String,char[]> seqlib = SeqLibrary.read(args[3]);
		ArrayList<String[]> pairs = PairReader.parse(args[4]);
		HashMap<String,CathScopEntry> cathscopinfo = CathScopHash.read(args[5]);
		
		
		double go = Double.parseDouble(args[0]);
		double ge = Double.parseDouble(args[1]);
//		double hbWeight;
//		double polWeight = Double.parseDouble(args[7]);
//		double ssWeight = Double.parseDouble(args[8]);
//		double substWeight = Double.parseDouble(args[9]);
		String file = args[6];
		double[][] substMatrix = QuasarMatrix.DAYHOFF_MATRIX;
		double[][] secStructMatrix = SecStructScores.matrix;
		double[][] polMatrix = scores.calcGotohInputMatrix(scores.calcPolarityScores());
		double[][] hbMatrix = scores.calcGotohInputMatrix(scores.calcHydropathyScores());
		
		
		AlignmentBenchmarker ab;
		BufferedWriter resultwriter;
		
		resultwriter = new BufferedWriter(new FileWriter(file));
		if(args[2].equals("freeshift")){
			gotoh = new FreeshiftMusterLite(go, ge,hbMatrix,polMatrix,secStructMatrix,substMatrix);
		}
		else if(args[2].equals("glocal")){
			gotoh = new GLocalMusterLite(go, ge, SeqLibrary.read(args[7]),hbMatrix,polMatrix,secStructMatrix,substMatrix,0.1,0.1,0.3,0.1);
		}
		else if(args[2].equals("global")){
			gotoh = new GlobalMusterLite(go, ge, SeqLibrary.read(args[7]),hbMatrix,polMatrix,secStructMatrix,substMatrix,0.1,0.1,0.3,0.1);
		}
	

		ab = new AlignmentBenchmarker(gotoh,seqlib,pairs,cathscopinfo,resultwriter);
		ab.benchmark();
		try {
			ab.printResults();
		} catch (IOException e) {
			System.out.println("cannot write output! (Benchmarker)");
		}
	}
}
