package bioinfo.alignment.kerbsch.temp;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.HashMap;
import java.util.Map.Entry;

public class DihedralAngles {

	public static void calc(String folderpath){
		HashMap<Character,HashMap<String,double[]>> angles = new HashMap<Character,HashMap<String,double[]>>();
		angles.put('A',new HashMap<String,double[]>());
		angles.put('C',new HashMap<String,double[]>());
		angles.put('D',new HashMap<String,double[]>());
		angles.put('E',new HashMap<String,double[]>());
		angles.put('F',new HashMap<String,double[]>());
		angles.put('G',new HashMap<String,double[]>());
		angles.put('H',new HashMap<String,double[]>());
		angles.put('I',new HashMap<String,double[]>());
		angles.put('K',new HashMap<String,double[]>());
		angles.put('L',new HashMap<String,double[]>());
		angles.put('M',new HashMap<String,double[]>());
		angles.put('N',new HashMap<String,double[]>());
		angles.put('P',new HashMap<String,double[]>());
		angles.put('Q',new HashMap<String,double[]>());
		angles.put('R',new HashMap<String,double[]>());
		angles.put('S',new HashMap<String,double[]>());
		angles.put('T',new HashMap<String,double[]>());
		angles.put('V',new HashMap<String,double[]>());
		angles.put('W',new HashMap<String,double[]>());
		angles.put('Y',new HashMap<String,double[]>());
		
		HashMap<Character,HashMap<String,Integer>> anglecounts = new HashMap<Character,HashMap<String,Integer>>();
		anglecounts.put('A',new HashMap<String,Integer>());
		anglecounts.put('C',new HashMap<String,Integer>());
		anglecounts.put('D',new HashMap<String,Integer>());
		anglecounts.put('E',new HashMap<String,Integer>());
		anglecounts.put('F',new HashMap<String,Integer>());
		anglecounts.put('G',new HashMap<String,Integer>());
		anglecounts.put('H',new HashMap<String,Integer>());
		anglecounts.put('I',new HashMap<String,Integer>());
		anglecounts.put('K',new HashMap<String,Integer>());
		anglecounts.put('L',new HashMap<String,Integer>());
		anglecounts.put('M',new HashMap<String,Integer>());
		anglecounts.put('N',new HashMap<String,Integer>());
		anglecounts.put('P',new HashMap<String,Integer>());
		anglecounts.put('Q',new HashMap<String,Integer>());
		anglecounts.put('R',new HashMap<String,Integer>());
		anglecounts.put('S',new HashMap<String,Integer>());
		anglecounts.put('T',new HashMap<String,Integer>());
		anglecounts.put('V',new HashMap<String,Integer>());
		anglecounts.put('W',new HashMap<String,Integer>());
		anglecounts.put('Y',new HashMap<String,Integer>());
		
		char precursor;
		char residue;
		char follower;
		double psi = 0.0;
		double phi = 0.0;
		String key = null;
				
		File dir = new File(folderpath);
		File[] files = dir.listFiles();
		int index = 0;
		for(File f : files){
			System.out.println(++index);
			try{
				BufferedReader in = new BufferedReader(new FileReader(f));
				String line;
				while((line = in.readLine()) != null){
					if(line.contains("NUMBER OF RESIDUES")){
						if(Integer.parseInt(line.split("\\s+")[1]) < 3){
							break;
						}
					}
					
					//read angles for every residue except of first and last residue
					//and store them in HashMap
					else if(line.trim().startsWith("#")){
						System.out.println(f.getName());
						
						//init steps
						line = in.readLine();						
						if(!line.matches("(.)*!(.)*") && line.charAt(13) != 'X'){	//chain border
							precursor = line.charAt(13);
							
							//lower-case letters in DSSP are Cysteines
							if(Character.isLowerCase(precursor)){
								precursor = 'C';
							}
							
						} else {
							precursor = 0;
						}
						line = in.readLine();
						if(!line.matches("(.)*!(.)*") && line.charAt(13) != 'X'){	//chain border
							residue = line.charAt(13);
							
							//lower-case letters in DSSP are Cysteines
							if(Character.isLowerCase(residue)){
								residue = 'C';
							}
							
							phi = Double.parseDouble(line.substring(103,109));
							psi = Double.parseDouble(line.substring(109,115));
						} else {
							residue = 0;
						}
						line = in.readLine();
						if(!line.matches("(.)*!(.)*") && line.charAt(13) != 'X'){	//chain border
							follower = line.charAt(13);
							
							//lower-case letters in DSSP are Cysteines
							if(Character.isLowerCase(follower)){
								follower = 'C';
							}
							
						} else {
							follower = 0;
						}
						if(precursor != 0 && residue != 0 && follower != 0){
							key = "" + precursor + follower;
							if(!angles.get(residue).containsKey(key)){
								angles.get(residue).put(key, new double[] {phi,psi});
								anglecounts.get(residue).put(key,new Integer(1));
							}else{
								angles.get(residue).get(key)[0] += phi;
								angles.get(residue).get(key)[1] += psi;
								anglecounts.get(residue).put(key, anglecounts.get(residue).get(key) + 1);
							}
							precursor = residue;
							residue = follower;
							phi = Double.parseDouble(line.substring(103,109));
							psi = Double.parseDouble(line.substring(109,115));
						}
						
						//read the rest of the DSSP
						while((line = in.readLine()) != null){							
							if(!line.matches("(.)*!(.)*") && line.charAt(13) != 'X'){	//chain border
								follower = line.charAt(13);
								
								//lower-case letters in DSSP are Cysteines
								if(Character.isLowerCase(follower)){
									follower = 'C';
								}
								
								if(precursor != 0 && residue != 0){
									key = "" + precursor + follower;
									if(!angles.get(residue).containsKey(key)){
										angles.get(residue).put(key, new double[] {phi,psi});
										anglecounts.get(residue).put(key,new Integer(1));
									}else{
										angles.get(residue).get(key)[0] += phi;
										angles.get(residue).get(key)[1] += psi;
										anglecounts.get(residue).put(key, anglecounts.get(residue).get(key) + 1);
									}
									precursor = residue;
									residue = follower;
									phi = Double.parseDouble(line.substring(103,109));
									psi = Double.parseDouble(line.substring(109,115));
								} else {
									precursor = residue;
									residue = follower;
								}
							}
						}
					}
				}
				in.close();
			} catch(IOException e){
				System.out.println("No secstruct input!");
			}
		}
		
		//write out results
		double phi_out;
		double psi_out;
		for (Entry<Character,HashMap<String,double[]>> amino: angles.entrySet()) {
//			System.out.println(amino.getKey()+": ______________");
			try{
				BufferedWriter out = new BufferedWriter(new FileWriter("angles/"+amino.getKey()+".angles"));
				for(Entry<String,double[]> pre_foll : amino.getValue().entrySet()){	
					phi_out = pre_foll.getValue()[0]/anglecounts.get(amino.getKey()).get(pre_foll.getKey());
					psi_out = pre_foll.getValue()[1]/anglecounts.get(amino.getKey()).get(pre_foll.getKey());
//					System.out.println(pre_foll.getKey() + ": " + pre_foll.getValue()[0] + " \\ " + anglecounts.get(amino.getKey()).get(pre_foll.getKey()));
//					System.out.println(pre_foll.getKey() + ": " + pre_foll.getValue()[1] + " \\ " + anglecounts.get(amino.getKey()).get(pre_foll.getKey()));
					out.append(pre_foll.getKey()+"\t"+phi_out+"\t"+psi_out+"\n");
				}
				out.close();
			} catch(IOException e){
				System.out.println("cannot write angles!");
			}
		}
	}
	
	public static HashMap<Character,HashMap<String,double[]>> read(String folderpath){
		BufferedReader in;
		HashMap<Character,HashMap<String,double[]>> angles = new HashMap<Character,HashMap<String,double[]>>();
		
		File dir = new File(folderpath);
		File[] files = dir.listFiles();
		
		char amino;
		String line = null;
		String[] temp;
		
		for(File f : files){
			try {
				in = new BufferedReader(new FileReader(f));
				amino = (f.getName()).charAt(0);
				angles.put(amino, new HashMap<String,double[]>());
				
				while((line = in.readLine()) != null){
					temp = line.split("\\s+");
					angles.get(amino).put(temp[0], new double[]{Double.parseDouble(temp[1]),Double.parseDouble(temp[2])});
				}
				in.close();
			} catch (IOException e) {
				System.out.println("cannot read angles");
			}
		}
		return angles;
	}
	
	public static void main(String[] args) {
		DihedralAngles.calc("DSSP");
	}

}
