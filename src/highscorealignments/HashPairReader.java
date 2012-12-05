package highscorealignments;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;

public class HashPairReader {
	
	public static HashMap<String,ArrayList<String>> readPairs(String infilename,String outfilename){
	HashMap<String,ArrayList<String>> pairs = new HashMap<String,ArrayList<String>>();

	try{
		BufferedReader in = new BufferedReader(new FileReader(infilename));
		String line;
		while((line = in.readLine()) != null){
			String[] temp = line.split("\\s");
			if(pairs.containsKey(temp[0])){
				pairs.get(temp[0]).add(temp[1]);
			} else{
				pairs.put(temp[0],new ArrayList<String>());
				pairs.get(temp[0]).add(temp[1]);	
			}
			if(pairs.containsKey(temp[1])){
				pairs.get(temp[1]).add(temp[0]);
			} else{
				pairs.put(temp[1],new ArrayList<String>());
				pairs.get(temp[1]).add(temp[0]);
			}
		}
		in.close();
		
		in = new BufferedReader(new FileReader(outfilename));
		while((line = in.readLine()) != null){
			String[] temp = line.split("\\s");
			if(pairs.containsKey(temp[0])){
				pairs.get(temp[0]).add(temp[1]);
			} else{
				pairs.put(temp[0],new ArrayList<String>());
				pairs.get(temp[0]).add(temp[1]);	
			}
			if(pairs.containsKey(temp[1])){
				pairs.get(temp[1]).add(temp[0]);
			} else{
				pairs.put(temp[1],new ArrayList<String>());
				pairs.get(temp[1]).add(temp[0]);
			}
		}
		
		in.close();
	} catch(IOException e){
		System.out.println("No pairs input!");
	}
	return pairs;
}
}
