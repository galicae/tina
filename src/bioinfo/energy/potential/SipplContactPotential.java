package bioinfo.energy.potential;

import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.List;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import bioinfo.proteins.AminoAcid;
import bioinfo.proteins.AminoAcidName;
import bioinfo.proteins.Atom;
import bioinfo.proteins.AtomType;
import bioinfo.proteins.PDBEntry;

public class SipplContactPotential extends AContactPotential{
	
	/**
	 * reads potential from VPOT file 
	 * e.g. vpot on biocluster
	 * @param filename
	 */
	public void readFromVPOTFile(String filename){
		BufferedReader br = null;
		try{
			int[] size = {6,6,26,26};
			initPotential(size);
			String line = null;
			Pattern firstlinePattern = Pattern.compile("#\\s+<k=(\\d+)\\s+(\\w{3})-(\\w{3})\\s+.*?>");
			Matcher firstlineMatcher = null;
			Pattern firstlinealterPattern = Pattern.compile("#\\s+<k=\\[(\\d+):\\d+\\]\\s+(\\w{3})-(\\w{3})\\s+.*?>");
			Matcher firstlinealterMatcher = null;
			Pattern contentlinePattern = Pattern.compile("(-?\\d+\\.\\d+)\\s+(-?\\d+\\.\\d+)\\s+(-?\\d+\\.\\d+)\\s+(-?\\d+\\.\\d+)\\s+(-?\\d+\\.\\d+)\\s+(-?\\d+\\.\\d+)");
			Matcher contentlineMatcher = null;
			br = new BufferedReader(new InputStreamReader(new FileInputStream(filename)));
			int k = 0;
			int i = 0;
			int j = 0;
			
			while((line = br.readLine()) != null){
				if(line.trim().startsWith("# <")){
					break;
				}
			}
			firstlineMatcher = firstlinePattern.matcher(line);
			firstlineMatcher.find();
			br.readLine();
			br.readLine();
			line = br.readLine();
			contentlineMatcher = contentlinePattern.matcher(line);
			contentlineMatcher.find();
			k = Integer.parseInt(firstlineMatcher.group(1))-2;
			i = AminoAcidName.getAAFromTLC(firstlineMatcher.group(2).trim()).getOneLetterCode().charAt(0)-65;
			j = AminoAcidName.getAAFromTLC(firstlineMatcher.group(3).trim()).getOneLetterCode().charAt(0)-65;
			for(int t = 0; t != 6; t++){
				int[] path = {k,t,i,j};
				potential.setValue(path,Double.parseDouble(contentlineMatcher.group(t+1)));
			}
			while((line = br.readLine()) != null){
				if(line.trim().startsWith("# <")){
					firstlineMatcher = firstlinePattern.matcher(line);
					if(firstlineMatcher.find()){
						br.readLine();
						br.readLine();
						line = br.readLine();
						contentlineMatcher = contentlinePattern.matcher(line);
						contentlineMatcher.find();
						k = Integer.parseInt(firstlineMatcher.group(1))-2;
						i = AminoAcidName.getAAFromTLC(firstlineMatcher.group(2).trim()).getOneLetterCode().charAt(0)-65;
						j = AminoAcidName.getAAFromTLC(firstlineMatcher.group(3).trim()).getOneLetterCode().charAt(0)-65;
						for(int t = 0; t != 6; t++){
							int[] path = {k,t,i,j};
							potential.setValue(path,Double.parseDouble(contentlineMatcher.group(t+1)));
						}
					} else {
						firstlinealterMatcher = firstlinealterPattern.matcher(line);
						if(firstlinealterMatcher.find()){
							br.readLine();
							br.readLine();
							line = br.readLine();
							contentlineMatcher = contentlinePattern.matcher(line);
							contentlineMatcher.find();
							k = Integer.parseInt(firstlinealterMatcher.group(1))-2;
							i = AminoAcidName.getAAFromTLC(firstlinealterMatcher.group(2).trim()).getOneLetterCode().charAt(0)-65;
							j = AminoAcidName.getAAFromTLC(firstlinealterMatcher.group(3).trim()).getOneLetterCode().charAt(0)-65;
							for(int t = 0; t != 6; t++){
								int[] path = {k,t,i,j};
								potential.setValue(path,Double.parseDouble(contentlineMatcher.group(t+1)));
							}
						}	
					}
				}
			}
		} catch (Exception e){
			e.printStackTrace();
		} finally {
			try {
				if (br != null) {
					br.close();
					br = null;
				}
			} catch (IOException e) {
				System.err.println("Error closing the vpot file: "+e.getLocalizedMessage());
				e.printStackTrace();
			}
		
		}
	}
	
	@Override
	public void calculateFromDATA(List<String> Ids, String dataLoc) {
		// nothing to do, due to the potential beeing precomputed
	}

	@Override
	public double scoreModel(PDBEntry model) {
		double energy = 0.0d;
		AminoAcid aa1 = null;
		Atom a1 = null;
		int p1 = 0;
		AminoAcid aa2 = null;
		Atom a2 = null;
		int p2 = 0;
		double dist = 0.0d;
		int d = 0;
		

		for(int l = 0; l != model.length()-7; l++){
			aa1 = model.getAminoAcid(l);
			a1 = aa1.getAtomByType(AtomType.CA);
			if(a1 == null){
				continue;
			}
			p1 = aa1.getName().getOneLetterCode().charAt(0)-65;
			for(int k = 0; k != 6; k++){
				aa2 = model.getAminoAcid(l+k+2);
				a2 = aa2.getAtomByType(AtomType.CA);
				if(a2 == null){
					continue;
				}
				p2 = aa2.getName().getOneLetterCode().charAt(0)-65;
				
				dist = Math.sqrt((a1.getPosition()[0]-a2.getPosition()[0])*(a1.getPosition()[0]-a2.getPosition()[0])+(a1.getPosition()[1]-a2.getPosition()[1])*(a1.getPosition()[1]-a2.getPosition()[1])+(a1.getPosition()[2]-a2.getPosition()[2])*(a1.getPosition()[2]-a2.getPosition()[2]));
								
				d = 0;
				if(dist >=3 && dist <= 15){
					for(int i = 5; i <= 13; i = i+2){
						if(dist <= i){
							break;
						}
						d++;
					}
					int[] path = {k,d,p1,p2};
					energy += potential.getByAddress(path).getValue();
				}
			}
		}
		return energy;
	}
	
	@Override
	public double[] getAminoScores(PDBEntry model){
		double[] scores = new double[model.length()];
		AminoAcid aa1 = null;
		Atom a1 = null;
		int p1 = 0;
		AminoAcid aa2 = null;
		Atom a2 = null;
		int p2 = 0;
		double dist = 0.0d;
		int d = 0;
		

		for(int l = 0; l != model.length()-7; l++){
			aa1 = model.getAminoAcid(l);
			a1 = aa1.getAtomByType(AtomType.CA);
			if(a1 == null){
				continue;
			}
			p1 = aa1.getName().getOneLetterCode().charAt(0)-65;
			for(int k = 0; k != 6; k++){
				aa2 = model.getAminoAcid(l+k+2);
				a2 = aa2.getAtomByType(AtomType.CA);
				if(a2 == null){
					continue;
				}
				p2 = aa2.getName().getOneLetterCode().charAt(0)-65;
				
				dist = Math.sqrt((a1.getPosition()[0]-a2.getPosition()[0])*(a1.getPosition()[0]-a2.getPosition()[0])+(a1.getPosition()[1]-a2.getPosition()[1])*(a1.getPosition()[1]-a2.getPosition()[1])+(a1.getPosition()[2]-a2.getPosition()[2])*(a1.getPosition()[2]-a2.getPosition()[2]));
								
				d = 0;
				if(dist >=3 && dist <= 15){
					for(int i = 5; i <= 13; i = i+2){
						if(dist <= i){
							break;
						}
						d++;
					}
					int[] path = {k,d,p1,p2};
					scores[l] = potential.getByAddress(path).getValue();
				}
			}
		}
		return scores;
	}
	
	@Override
	public double getNativeScoring(){
		return scoreModel(this.model);
	}
	
	@Override
	public double getSequenceScoring(AminoAcidName[] sequence){
		AminoAcid[] aas = new AminoAcid[this.model.length()];
		for(int i = 0; i != this.model.length(); i++){
			aas[i] = new AminoAcid(sequence[i],i, model.getAminoAcid(i).getAtoms());
		}
		PDBEntry newModel = new PDBEntry(model.getId()+model.getChainID()+String.format("%02d",model.getChainIDNum())+"_shuffled", aas);
		return scoreModel(newModel);
	}
	
	@Override
	public double[] getNativeAminoScoring(){
		return getAminoScores(this.model);
	}
	
	@Override
	public double[] getSequenceAminoScoring(AminoAcidName[] sequence){
		AminoAcid[] aas = new AminoAcid[this.model.length()];
		for(int i = 0; i != this.model.length(); i++){
			aas[i] = new AminoAcid(sequence[i],i, model.getAminoAcid(i).getAtoms());
		}
		PDBEntry newModel = new PDBEntry(model.getId()+model.getChainID()+String.format("%02d",model.getChainIDNum())+"_shuffled", aas);
		return getAminoScores(newModel);
	}	
	

	public static void main(String[] args) {
		SipplContactPotential pot = new SipplContactPotential();
		//pot.readFromVPOTFile("/Users/andreseitz/Desktop/sippl.vpot");
		//pot.writeToFile("/Users/andreseitz/Desktop/sippl.pot");
		pot.readFromFile("/Users/andreseitz/Desktop/sippl.pot");
		
		
		System.out.println("did magic");
		
	}
}
