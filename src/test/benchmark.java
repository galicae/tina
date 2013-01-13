package test;

import java.io.IOException;
import java.util.HashMap;

import bioinfo.alignment.kerbsch.AlignmentBenchmarker;
import bioinfo.alignment.kerbsch.SecStructReader;

public class benchmark {
	public static void main(String[] args) {
		HashMap<String,char[]> seqlib = SecStructReader.read("../Gobi_old/DSSP");
		AlignmentBenchmarker ab = new AlignmentBenchmarker(args,seqlib);
		ab.benchmark();
		try {
			ab.printResults();
		} catch (IOException e) {
			System.out.println("cannot write output! (Benchmarker)");
		}
	}
}
