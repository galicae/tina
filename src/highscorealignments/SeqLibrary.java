package highscorealignments;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.HashMap;

public class SeqLibrary {
			
	public static HashMap<String, char[]> parse(String filename){
		HashMap<String, char[]> sequences = new HashMap<String, char[]>();
		try{
			BufferedReader in = new BufferedReader(new FileReader(filename));
			String line;
			while((line = in.readLine()) != null){
				String[] temp = line.split(":");
				sequences.put(temp[0], temp[1].toCharArray());
			}
			in.close();
		} catch(IOException e){
			System.out.println("No seqlib input!");
		}
		return sequences;
	}
}

