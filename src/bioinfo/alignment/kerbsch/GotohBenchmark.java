package bioinfo.alignment.kerbsch;

import highscorealignments.CathScopEntry;
import highscorealignments.CathScopHash;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;

import bioinfo.alignment.Aligner;
import bioinfo.alignment.gotoh.FreeshiftSequenceGotoh;
import bioinfo.alignment.gotoh.GlobalSequenceGotoh;
import bioinfo.alignment.gotoh.LocalSequenceGotoh;
import bioinfo.alignment.kerbsch.temp.InitClass;
import bioinfo.alignment.kerbsch.temp.PairReader;
import bioinfo.alignment.kerbsch.temp.SecStructScores;
import bioinfo.alignment.kerbsch.temp.SeqLibrary;
import bioinfo.alignment.matrices.QuasarMatrix;

public class GotohBenchmark {

	public static void main(String[] args) throws IOException {
		InitClass scores = new InitClass();
		Aligner gotoh = null;
		BufferedWriter resultwriter;

//		double go = Double.parseDouble(args[0]);
//		double ge = Double.parseDouble(args[1]);
		
		double[][] matrix;
		HashMap<String,char[]> seqlib = SeqLibrary.read(args[4]);
		
		if(args[2].equals("dayhoff")){
			matrix = QuasarMatrix.DAYHOFF_MATRIX;
		} else if(args[2].equals("polarity")){
			matrix = scores.calcGotohInputMatrix(scores.calcPolarityScores());
		} else if(args[2].equals("hydrophob")){
			matrix = scores.calcGotohInputMatrix(scores.calcHydropathyScores());
		} else {
			matrix = SecStructScores.matrix;
			seqlib = SeqLibrary.read(args[4]);
		}
		

		
		ArrayList<String[]> pairs = PairReader.parse(args[5]);
		HashMap<String,CathScopEntry> cathscopinfo = CathScopHash.read(args[6]);
		
		
		for (int go = -8; go > -16; go--) {
			for (int ge = -1; ge >= go ; ge--) {
				if(args[3].equals("freeshift")){
					gotoh = new FreeshiftSequenceGotoh(go, ge, matrix);
				}
				else if(args[3].equals("local")){
					gotoh = new LocalSequenceGotoh(go, ge, matrix);
				}
				else if(args[3].equals("global")){
					gotoh = new GlobalSequenceGotoh(go, ge, matrix);
				}
				else if(args[3].equals("glocal")){
					gotoh = new GLocalSequenceGotoh(go, ge, matrix);
				}
				resultwriter = new BufferedWriter(new FileWriter("../dayhoff/"+args[2]+"_"+args[3]+"_"+go+"_"+ge+".bm"));
				AlignmentBenchmarker ab = new AlignmentBenchmarker(gotoh,seqlib,pairs,cathscopinfo,resultwriter);
				ab.benchmark();
				try {
					ab.printResults();
				} catch (IOException e) {
					System.out.println("cannot write output! (Benchmarker)");
				}
			}
		}
	}

}
