package highscorealignments;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;

public class HashPairReaderCathScop {
	
	public static HashMap<String,ArrayList<CathScopEntry[]>> readPairs(String filename){
	HashMap<String,ArrayList<CathScopEntry[]>> pairs = new HashMap<String,ArrayList<CathScopEntry[]>>();

	try{
		BufferedReader in = new BufferedReader(new FileReader(filename));
		String line;
		String lastid = "";
		String actid;
		while((line = in.readLine()) != null){
			String[] temp = line.split("\\s");
			actid = temp[0];
			if(actid.equals(lastid)){
				pairs.get(actid).add(new CathScopEntry[2]);
				String[] cath1 = temp[temp.length-2].split("\\.");
				String[] cath2 = temp[temp.length-1].split("\\.");
				String[] scop1 = temp[temp.length-6].split("\\.");
				String[] scop2 = temp[temp.length-5].split("\\.");
				pairs.get(actid).get(pairs.get(actid).size()-1)[0] = new CathScopEntry(temp[0],cath1[0].charAt(0),Integer.parseInt(cath1[3]),Integer.parseInt(cath1[2]),Integer.parseInt(cath1[1]),
						Integer.parseInt(scop1[0]),Integer.parseInt(scop1[3]),Integer.parseInt(scop1[2]),Integer.parseInt(scop1[1]));
				pairs.get(actid).get(pairs.get(actid).size()-1)[1] = new CathScopEntry(temp[1],cath2[0].charAt(0),Integer.parseInt(cath2[3]),Integer.parseInt(cath2[2]),Integer.parseInt(cath2[1]),
						Integer.parseInt(scop2[0]),Integer.parseInt(scop2[3]),Integer.parseInt(scop2[2]),Integer.parseInt(scop2[1]));		
			} else{
				lastid = actid;
				pairs.put(actid,new ArrayList<CathScopEntry[]>());
				pairs.get(actid).add(new CathScopEntry[2]);
				String[] cath1 = temp[temp.length-2].split("\\.");
				String[] cath2 = temp[temp.length-1].split("\\.");
				String[] scop1 = temp[temp.length-6].split("\\.");
				String[] scop2 = temp[temp.length-5].split("\\.");
				pairs.get(actid).get(pairs.get(actid).size()-1)[0] = new CathScopEntry(temp[0],cath1[0].charAt(0),Integer.parseInt(cath1[3]),Integer.parseInt(cath1[2]),Integer.parseInt(cath1[1]),
						Integer.parseInt(scop1[0]),Integer.parseInt(scop1[3]),Integer.parseInt(scop1[2]),Integer.parseInt(scop1[1]));
				pairs.get(actid).get(pairs.get(actid).size()-1)[1] = new CathScopEntry(temp[1],cath2[0].charAt(0),Integer.parseInt(cath2[3]),Integer.parseInt(cath2[2]),Integer.parseInt(cath2[1]),
						Integer.parseInt(scop2[0]),Integer.parseInt(scop2[3]),Integer.parseInt(scop2[2]),Integer.parseInt(scop2[1]));
		
			}
		}
		in.close();
	} catch(IOException e){
		System.out.println("No pairs input!");
	}
	return pairs;
}
}
