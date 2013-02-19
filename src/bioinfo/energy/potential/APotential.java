package bioinfo.energy.potential;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.InputStreamReader;
import java.io.OutputStreamWriter;
import java.util.HashMap;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

public abstract class APotential implements IEnergy{
	
	PotentialDimension<Double> potential;
	
	@Override
	public void initPotential(int[] size) {
		potential = new PotentialDimension<Double>(size, 0.0);
	}
	
	@Override
	public void writeToFile(String filename) {
		HashMap<Integer,Integer> size = potential.getSize();
		BufferedWriter bw;
		try {
			bw = new BufferedWriter(new OutputStreamWriter(new FileOutputStream(filename)));
			bw.append("<header>");
			bw.append("\n");
			bw.append("<levels count=\""+size.size()+"\">");
			bw.append("\n");
			for(int l = 0; l != size.size(); l++){
				if(size.containsKey(l)){
					bw.append("<level number=\""+l+"\" size=\""+size.get(l)+"\"/>");
				}else{
					bw.append("<level number=\""+l+"\" size=\""+1+"\"/>");
				}
				bw.append("\n");
			}
			bw.append("</levels>");
			bw.append("\n");
			bw.append("</header>");
			bw.append("\n");
			
			potential.writeAsXML(bw);
			bw.flush();
			bw.close();
		} catch(Exception e){
			e.printStackTrace();
		}
	}
	
	@SuppressWarnings("unused")
	@Override
	public void readFromFile(String filename) {
		try{
			BufferedReader br = new BufferedReader(new InputStreamReader(new FileInputStream(filename)));
			String line = null;
			
			Pattern headerPattern = Pattern.compile("<header>");
			Matcher headerMatcher = null;
			Pattern levelsPattern = Pattern.compile("<levels count=\"(\\d+)\">");
			Matcher levelsMatcher = null;
			Pattern levelPattern = Pattern.compile("<level number=\"(\\d+)\" size=\"(\\d+)\"/>");
			Matcher levelMatcher = null;
			Pattern levelsClosePattern = Pattern.compile("</levels>");
			Matcher levelsCloseMatcher = null;
			Pattern headerClosePattern = Pattern.compile("</header>");
			Matcher headerCloseMatcher = null;
			
			
			//int levels = 0; // needed for proofreading
			HashMap<Integer,Integer> size = new HashMap<Integer,Integer>();
			
			while((line = br.readLine()) != null){
				if((headerMatcher = headerPattern.matcher(line)).find()){
					while((line = br.readLine()) != null){
						if((levelsMatcher = levelsPattern.matcher(line)).find()){
							//levels = Integer.parseInt(levelsMatcher.group(1));
							while((line = br.readLine()) != null){
								if((levelMatcher = levelPattern.matcher(line)).find()){
									size.put(Integer.parseInt(levelMatcher.group(1)),Integer.parseInt(levelMatcher.group(2)));
								}else if((levelsCloseMatcher = levelsClosePattern.matcher(line)).find()){
									break;
								}
							}
							if((levelsCloseMatcher = levelsClosePattern.matcher(line)).find()){
								break;
							}
						}else if((levelsCloseMatcher = levelsClosePattern.matcher(line)).find()){
							break;
						}
					}
					if((headerCloseMatcher = headerClosePattern.matcher(line)).find()){
						break;
					}
				}else if((headerCloseMatcher = headerClosePattern.matcher(line)).find()){
					break;
				}
			}
			
			
			//proofreading
			//if(levels == size.size()){...}
			int[] sizes = new int[size.size()];
			int[] pointer = new int[size.size()];
			for(int i = 0; i != size.size(); i++){
				pointer[i] = -1;
				sizes[i] = size.get(i);
			}
			
			
			potential = new PotentialDimension<Double>(sizes,0.0);
			
			Pattern dimensionPattern = Pattern.compile("<dimension level=\"(\\d+)\">");
			Matcher dimensionMatcher = null;
			Pattern dimensionClosePattern = Pattern.compile("</dimension>");
			Matcher dimensionCloseMatcher = null;
			Pattern dataPattern = Pattern.compile("<data value=\"(.*?)\"/>");
			Matcher dataMatcher = null;
			int levelCt = -1;
			
			while((line = br.readLine()) != null){
				if((dimensionMatcher = dimensionPattern.matcher(line)).find()){
					levelCt++;
					//proofreading
					//if(Integer.parseInt(dimensionMatcher.group(1)) == levelCt){...}
					if(levelCt >= 0){
						pointer[levelCt]++;
					}
				}else if((dimensionCloseMatcher = dimensionClosePattern.matcher(line)).find()){
					pointer[levelCt] = -1;
					levelCt--;
					if(levelCt >= 0){
						pointer[levelCt]++;
					}else{
						break;
					}
				}else if((dataMatcher = dataPattern.matcher(line)).find()){
					potential.setValue(pointer,Double.parseDouble(dataMatcher.group(1)));
					pointer[levelCt]++;
				}
			}
			
				
			br.close();
			
			
		}catch(Exception e){
			e.printStackTrace();
		}
	}

}
