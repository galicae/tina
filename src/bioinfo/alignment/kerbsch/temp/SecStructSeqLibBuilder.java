package bioinfo.alignment.kerbsch.temp;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.HashMap;
import java.util.Map.Entry;


public class SecStructSeqLibBuilder {
	public static void build(HashMap<String,char[]> seqlib, String outpath){
		BufferedWriter out;
//		HashMap<String,char[]> seqlib_ref = SeqLibrary.read("../GoBi/referenz/domains.seqlib");
		try{
			out = new BufferedWriter(new FileWriter(outpath));
			for(Entry<String,char[]> entry : seqlib.entrySet()){
//				if(seqlib_ref.containsKey(entry.getKey())){
				out.write(entry.getKey() + ":");
				out.write(entry.getValue());
				out.write("\n");
//				}
			}
			out.close();
		}catch(IOException e){
			System.out.println("cannot write file(SecStructBuilder)");
		}
	}
	
//	public static void main (String[] args){
//		HashMap<String,char[]> seqlib = SSReaderFromDir.read("../GoBi/DSSP");
//		SecStructSeqLibBuilder.build(seqlib, "secstruct.seqlib");
//	}
}
