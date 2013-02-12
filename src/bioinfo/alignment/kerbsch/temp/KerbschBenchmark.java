package bioinfo.alignment.kerbsch.temp;

import highscorealignments.CathScopEntry;
import highscorealignments.CathScopHash;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;

import bioinfo.alignment.kerbsch.AlignmentBenchmarker;
import bioinfo.alignment.kerbsch.Kerbsch;
import bioinfo.alignment.matrices.QuasarMatrix;

public class KerbschBenchmark {


	public static void main(String[] args) throws IOException {
		InitClass init = new InitClass();
		double[][] seqMatrix = QuasarMatrix.DAYHOFF_MATRIX;
		double[][] polMatrix = init.calcGotohInputMatrix(init.calcPolarityScores());
		double[][] hbMatrix = init.calcGotohInputMatrix(init.calcHydropathyScores());
		double[][] secStructMatrix = SecStructScores.matrix;
		
		BufferedWriter out = new BufferedWriter(new FileWriter("kerbsch.bench"));
		HashMap<String,char[]> seclib = SeqLibrary.read("../tina/secstruct.seqlib");
		HashMap<String,char[]> seqlib = SeqLibrary.read("../tina/domains.seqlib");
		ArrayList<String[]> pairs = PairReader.parse("../tina/410List.pairs");
		HashMap<String,CathScopEntry> cathscopinfo = CathScopHash.read("../tina/cathscop.seqlib");
		
		Kerbsch kerbsch = new Kerbsch(200.0,50.0,seqMatrix,hbMatrix,polMatrix,secStructMatrix,seclib);
		AlignmentBenchmarker bench = new AlignmentBenchmarker(kerbsch,seqlib,pairs,cathscopinfo,out);
		
		bench.benchmark();
		bench.printResults();
	}

}
