package bioinfo.proteins;

import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStreamReader;

public class SSCCReader {
	
	public static SSCCEntry readSSCC(String filepath){
		BufferedReader in;
		String line;
		String id = null;
		String[] temp;
		SSCCLine[] sscclines = null;
		int index = -1;
		
		try{
			in = new BufferedReader(new InputStreamReader(new FileInputStream(filepath)));
			
			//initialize
			temp = in.readLine().split("\\s+");
			id = temp[0];
			sscclines = new SSCCLine[Integer.parseInt(temp[1])];
			
			//readin
			while((line = in.readLine()) != null){
				temp = line.split("\\s+");
				index++;
				sscclines[index] = new SSCCLine(AminoAcidName.getAAFromOLC(temp[0]),temp[1].charAt(0),Integer.parseInt(temp[2]),Integer.parseInt(temp[3]));
			}
			
			in.close();
		} catch(IOException e){
			System.out.println("cannot read from sscc-file");
		}
		
		//return
		SSCCEntry entry = new SSCCEntry(id,sscclines);
		return entry;
	}
}
