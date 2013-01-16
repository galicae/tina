package bioinfo.alignment.kerbsch;

import java.io.IOException;
import java.util.HashMap;


public class benchmark {
	public static void main(String[] args) {
		HashMap<String,char[]> seqlib = SeqLibrary.read(args[3]);
		AlignmentBenchmarker ab = new AlignmentBenchmarker(args,seqlib);
		ab.benchmark();
//		try {
//			ab.printResults();
//		} catch (IOException e) {
//			System.out.println("cannot write output! (Benchmarker)");
//		}
	}
}
