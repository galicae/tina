package bioinfo.alignment.kerbsch.temp;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.HashMap;
import java.util.Map.Entry;

import bioinfo.proteins.DSSPEntry;
import bioinfo.proteins.SecStructEight;

import db.mysql.DBConnector;
import db.mysql.LocalConnection;


public class SecStructSeqLibBuilder {
	public static void main(String[] args){
		LocalConnection con = new LocalConnection();
		DBConnector dbcon = new DBConnector(con);	
		BufferedWriter out;
		HashMap<String,char[]> seqlib = SeqLibrary.read("../full_domains.seqlib");
		DSSPEntry dsspentry;
		
		try{
			out = new BufferedWriter(new FileWriter("../test_secstruct.seqlib"));
			
			for(Entry<String,char[]> entry : seqlib.entrySet()){
				dsspentry = dbcon.getDSSP(entry.getKey());
				out.write(entry.getKey() + ":");
				for(SecStructEight ss : dsspentry.getSecondaryStructure()){
					out.append(ss.getThreeClassAnalogon().getCharRepres());
				}
				out.write("\n");
			}
			out.close();
		}catch(IOException e){
			System.out.println("cannot write file(SecStructBuilder)");
		}
	}
}
