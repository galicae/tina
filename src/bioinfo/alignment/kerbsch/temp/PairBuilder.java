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
import java.util.List;
import java.util.Map.Entry;

public class PairBuilder {
	
	public static void buildSmallSubsetPairs(){
		BufferedReader in;
		BufferedWriter out;
		ArrayList<String> ids = new ArrayList<String>();
		
		try{
			in = new BufferedReader(new FileReader("../410.list"));
			String line;
			while((line = in.readLine()) != null){
				ids.add(line.trim());
			}
			in.close();
			
			List<String> pairs = new ArrayList<String>();
			
			out = new BufferedWriter(new FileWriter("410List.unique_pairs"));
			for (int i = 0; i < ids.size(); i++) {
				for (int j = 0; j < ids.size(); j++) {
					pairs.add(ids.get(i)+ids.get(j));
					if(ids.get(i) != ids.get(j) && !pairs.contains(ids.get(j)+ids.get(i))){		
						out.append(ids.get(i) + "\t" + ids.get(j) + "\n");
					}
				}
			}
			out.close();
		} catch (IOException e) {
			System.out.println("cannot read voretargets");
		}
	}
	
	
	public static void buildVorePairs(){
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
		PairBuilder.buildSmallSubsetPairs();
	}
}
