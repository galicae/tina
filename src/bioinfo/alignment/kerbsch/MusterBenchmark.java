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
		BufferedWriter resultwriter = new BufferedWriter(new FileWriter(args[6]));
		
		double go = Double.parseDouble(args[0]);
		double ge = Double.parseDouble(args[1]);
		double hbWeight = Double.parseDouble(args[7]);
		double polWeight = Double.parseDouble(args[8]);
		double ssWeight = Double.parseDouble(args[9]);
		double substWeight = Double.parseDouble(args[10]);
		double[][] substMatrix = QuasarMatrix.DAYHOFF_MATRIX;
		double[][] secStructMatrix = SecStructScores.matrix;
		double[][] polMatrix = scores.calcGotohInputMatrix(scores.calcPolarityScores());
		double[][] hbMatrix = scores.calcGotohInputMatrix(scores.calcHydropathyScores());
		
		if(args[2].equals("freeshift")){
			gotoh = new FreeshiftMusterLite(go, ge,hbMatrix,polMatrix,secStructMatrix,substMatrix);
		}
		else if(args[2].equals("local")){
//			gotoh = new LocalSequenceGotoh(go, ge,hbMatrix,polMatrix,secStructMatrix,substMatrix);
		}
		else if(args[2].equals("global")){
			gotoh = new GlobalMusterLite(go, ge, SeqLibrary.read(args[11]),hbMatrix,polMatrix,secStructMatrix,substMatrix,hbWeight,polWeight,ssWeight,substWeight);
		}
		
		HashMap<String,char[]> seqlib = SeqLibrary.read(args[3]);
		ArrayList<String[]> pairs = PairReader.parse(args[4]);
		HashMap<String,CathScopEntry> cathscopinfo = CathScopHash.read(args[5]);
	
		AlignmentBenchmarker ab = new AlignmentBenchmarker(gotoh,seqlib,pairs,cathscopinfo,resultwriter);
		ab.benchmark();
		try {
			ab.printResults();
		} catch (IOException e) {
			System.out.println("cannot write output! (Benchmarker)");
		}
	}
}
