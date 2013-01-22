package bioinfo.alignment.kerbsch.temp;

import highscorealignments.CathScopEntry;
import highscorealignments.CathScopHash;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map.Entry;

public class VoreScorePairBuilder {
	public static void build(){
		BufferedReader in;
		BufferedWriter out;
		HashMap<String,CathScopEntry> cathscopinfo = CathScopHash.read("cathscop.seqlib");
		HashMap<String,char[]> seqlib = SeqLibrary.read("../GoBi/referenz/domains.seqlib");
		ArrayList<String> targets = new ArrayList<String>();
		
		try{
			in = new BufferedReader(new FileReader("vorescore.foldtest.targets"));
			String line;
			while((line = in.readLine()) != null){
				targets.add(line.split("\\s+")[0]);
			}
			in.close();
		} catch (IOException e) {
			System.out.println("cannot read voretargets");
		}
		
		int index = 0;
		try {
			System.out.println(seqlib.size());
			out = new BufferedWriter(new FileWriter("vorescore.pairs"));
			for(String s : targets){
				for(Entry<String,char[]> entry : seqlib.entrySet()){
					if(!s.equals(entry.getKey())){
						out.append(s + "\t" + entry.getKey() + "\n");
					}
				}
			}
			out.close();
		} catch (IOException e) {
			System.out.println("cannot write vorescorepairs");
		}
	}
	
	public static void main(String[] args){
		VoreScorePairBuilder.build();
	}
}
