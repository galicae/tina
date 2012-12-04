package highscorealignments;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;

public class PairReaderCathScop {
	
	public static ArrayList<CathScopEntry[]> readPairs(String filename){
	ArrayList<CathScopEntry[]> pairs = new ArrayList<CathScopEntry[]>();

		try{
			BufferedReader in = new BufferedReader(new FileReader(filename));
			String line;
			while((line = in.readLine()) != null){
				String[] temp = line.split("\\s");
				String[] cath1entry = temp[temp.length-2].split("\\.");
				String[] cath2entry = temp[temp.length-1].split("\\.");
				String[] scop1entry = temp[temp.length-6].split("\\.");
				String[] scop2entry = temp[temp.length-5].split("\\.");
				pairs.add(new CathScopEntry[2]);
				pairs.get(pairs.size()-1)[0] = new CathScopEntry(temp[0],cath1entry[0].charAt(0),Integer.parseInt(cath1entry[3]),Integer.parseInt(cath1entry[2]),Integer.parseInt(cath1entry[1]),
						Integer.parseInt(scop1entry[0]),Integer.parseInt(scop1entry[3]),Integer.parseInt(scop1entry[2]),Integer.parseInt(scop1entry[1]));
				pairs.get(pairs.size()-1)[1] = new CathScopEntry(temp[1],cath2entry[0].charAt(0),Integer.parseInt(cath2entry[3]),Integer.parseInt(cath2entry[2]),Integer.parseInt(cath2entry[1]),
						Integer.parseInt(scop2entry[0]),Integer.parseInt(scop2entry[3]),Integer.parseInt(scop2entry[2]),Integer.parseInt(scop2entry[1]));
			}
			in.close();
		} catch(IOException e){
			System.out.println("No pairs input!");
		}
		return pairs;
	}
}
