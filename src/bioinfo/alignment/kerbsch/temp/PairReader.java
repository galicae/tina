package bioinfo.alignment.kerbsch.temp;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;

public class PairReader {
		
	public static ArrayList<String[]> parse(String filename){
		ArrayList<String[]> pairs = new ArrayList<String[]>();
		try{
			BufferedReader in = new BufferedReader(new FileReader(filename));
			String line;
			while((line = in.readLine()) != null){
				String[] temp = line.split("\\s");
				pairs.add(new String[]{temp[0],temp[1]});
			}
			in.close();
		} catch(IOException e){
			System.out.println("No pairs input!");
		}
		return pairs;
	}
}

