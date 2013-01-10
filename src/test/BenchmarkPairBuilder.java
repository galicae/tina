package test;

import highscorealignments.HashPairReader;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map.Entry;

public class BenchmarkPairBuilder {

	public static void main(String[] args) throws IOException {
		BufferedWriter out;
		HashMap<String,ArrayList<String>> pairs = HashPairReader.readPairs("../Gobi_old/referenz/cathscop.inpairs", "../Gobi_old/referenz/cathscop.outpairs");
		
		out = new BufferedWriter(new FileWriter("benchmark.pairs"));
		for(Entry<String,ArrayList<String>> e : pairs.entrySet()){
			for(String s : e.getValue()){
				out.append(e.getKey()+"\t"+s+"\n");
			}
		}
		out.close();
		System.out.println("done!");
	}

}
