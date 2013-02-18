package bioinfo.proteins.fragm3nt.run;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.util.HashMap;
import java.util.LinkedList;

import static bioinfo.proteins.fragm3nt.run.RunHelper.*;

/**
 * purpose of this class is to evaluate the TM score of the models produced by
 * the AlignmentAssembler class on the bucket05 data
 * 
 * @author galicae
 * 
 */
public class AliAssVsGotoh {
	static int fragLength = 9;
	static int extent = 4;
	static String desktop = "/home/galicae/Desktop/STRUCTURES/";
	static BufferedWriter resultWriter;

	public static void main(String[] args) throws Exception {
		BufferedWriter resultWriter  = new BufferedWriter(new FileWriter("spnf_gotohVSass"));
		// bucket file direction
		String bucket = "spnf";// args[0];
		// first find how many different sequences there are in the bucket
		// so first read bucket file
		// idea: use every ID as key and map to it a list with all the sequence
		// IDs that are aligned against the first ID
		HashMap<String, LinkedList<String>> map = new HashMap<String, LinkedList<String>>();
		BufferedReader r = new BufferedReader(new FileReader(bucket));

		String line = "";
		String[] lineArr = new String[2];
		while ((line = r.readLine()) != null) {
			lineArr = line.split(" ");
			if (!map.containsKey(lineArr[0])) {
				map.put(lineArr[0], new LinkedList<String>());
				map.get(lineArr[0]).add(lineArr[0]);
				map.get(lineArr[0]).add(lineArr[1]);
			} else {
				map.get(lineArr[0]).add(lineArr[1]);
			}
			if (!map.containsKey(lineArr[1])) {
				map.put(lineArr[1], new LinkedList<String>());
				map.get(lineArr[1]).add(lineArr[1]);
				map.get(lineArr[1]).add(lineArr[0]);
			} else {
				map.get(lineArr[1]).add(lineArr[0]);
			}
		}
		r.close();

		// now for every entry do the stuff from AlignmentAssemblyTest
		for (String s : map.keySet()) {
			System.out.println(s);
			try {
				// calculate Gotoh score
				resultWriter.write(s + "\t");
				
				String seq1ID = map.get(s).get(0);
				String seq2ID = map.get(s).get(1);
				double gotScore = gotohPart(seq1ID, seq2ID);
				System.out.print(gotScore + "\t");
				resultWriter.write(gotScore + "\t");
				// calculate Fragm3nt score
				
				double assScore = doMagic(map.get(s));
				System.out.println(assScore);
				resultWriter.write(assScore + "\n");
			} catch (Exception e) {
				e.printStackTrace();
				resultWriter.write("#\n");
//				System.exit(0);
				continue;
			}
		}
		resultWriter.close();
	}	
}
