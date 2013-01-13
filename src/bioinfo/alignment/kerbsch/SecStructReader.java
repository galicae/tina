package bioinfo.alignment.kerbsch;


import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.HashMap;

public class SecStructReader {
	
	public static HashMap<String,char[]> read(String folderpath){
		HashMap<String,char[]> secstruct = new HashMap<String,char[]>();
		char residue;
		char[] struct = null;	
		int rescount;
		int index;
		File dir = new File(folderpath);
		File[] files = dir.listFiles();
		
		for(File f : files){
			String id = f.getName().split("\\.")[0];
			index = 0;
			try{
				BufferedReader in = new BufferedReader(new FileReader(f));
				String line;
				while((line = in.readLine()) != null){
					
					//find number of residues
					if(line.contains("TOTAL NUMBER OF RESIDUES")){
						rescount = Integer.parseInt(line.split("\\s+")[1]);
						struct = new char[rescount];
					}
					
					//read secondary structure
					else if(line.trim().contains("#")){
						while((line = in.readLine()) != null){							
							if(!line.matches("(.)*!(.)*")){	//chain border
								residue = line.charAt(16);
								if(residue == 'H' || residue == 'G' || residue == 'I'){
									struct[index] = 'H';
								}
								else if(residue == 'B' || residue == 'E'){
									struct[index] = 'E';
	
								}
								else{
									struct[index] = 'C';
								}
								index++;
							}
						}
					}
				}
				secstruct.put(id,struct);
				in.close();
			} catch(IOException e){
				System.out.println("No secstruct input!");
			}
		}
		return secstruct;
	}
	
}
