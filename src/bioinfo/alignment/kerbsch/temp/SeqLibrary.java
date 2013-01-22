package bioinfo.alignment.kerbsch.temp;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.HashMap;

public class SeqLibrary {
			
	public static HashMap<String, char[]> read(String filename){
		HashMap<String, char[]> sequences = new HashMap<String, char[]>();
		try{
			BufferedReader in = new BufferedReader(new FileReader(filename));
			String line;
			String[] temp;
			while((line = in.readLine()) != null){
				temp = line.split(":");
				sequences.put(temp[0], temp[1].toCharArray());
			}
			in.close();
		} catch(IOException e){
			System.out.println("No seqlib input!");
		}
		return sequences;
	}
}

