package bioinfo.alignment.kerbsch.temp;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;

import bioinfo.alignment.kerbsch.Kerbsch;
import bioinfo.alignment.kerbsch.ValidateAlignments;
import bioinfo.alignment.matrices.QuasarMatrix;

public class KerbschBenchmark {


	public static void main(String[] args) throws IOException {
		InitClass init = new InitClass();
		double[][] seqMatrix = QuasarMatrix.DAYHOFF_MATRIX;
		double[][] polMatrix = init.calcGotohInputMatrix(init.calcPolarityScores());
		double[][] hbMatrix = init.calcGotohInputMatrix(init.calcHydropathyScores());
		double[][] secStructMatrix = SecStructScores.matrix;
		
		BufferedWriter out = new BufferedWriter(new FileWriter("kerbsch.bench"));
		HashMap<String,char[]> seclib = SeqLibrary.read("../full_secstruct.seqlib");
		HashMap<String,char[]> seqlib = SeqLibrary.read("../full_domains.seqlib");
		
		String targetsfolder = "../QUERY2TEMPLATELIST";
		BufferedReader in = new BufferedReader(new FileReader("../small.list"));
		String line;
		ArrayList<String> targetlist = new ArrayList<String>();
		while((line = in.readLine())!= null){
			targetlist.add(line);
		}
		in.close();
		
		Kerbsch kerbsch = new Kerbsch(-10.0,-1.0,seqMatrix,hbMatrix,polMatrix,secStructMatrix,seclib);
		ValidateAlignments bench = new ValidateAlignments(kerbsch,seqlib,targetlist,targetsfolder,out);
		
		bench.benchmark();
		bench.printResults();
	}

}
