package bioinfo.alignment.kerbsch;

import highscorealignments.CathScopEntry;
import highscorealignments.CathScopHash;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
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
		

		
		BufferedReader in = new BufferedReader(new FileReader(args[5]));
		String line;
		ArrayList<String> targets = new ArrayList<String>();
		while((line = in.readLine())!= null){
			targets.add(line);
		}
		in.close();
		String targetsfolder = args[6];
		
		
		for (int go = -1; go > -6; go--) {
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
				resultwriter = new BufferedWriter(new FileWriter("../"+args[2]+"/"+args[2]+"_"+go+"_"+ge+".bm"));
				ValidateAlignments va = new ValidateAlignments(gotoh,seqlib,targets,targetsfolder,resultwriter);
				va.benchmark();
				va.printResults();
			}
		}
	}

}
