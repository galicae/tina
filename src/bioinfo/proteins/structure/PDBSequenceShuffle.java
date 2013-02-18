package bioinfo.proteins.structure;

import java.util.ArrayList;
import java.util.List;
import java.util.Random;

import bioinfo.proteins.AminoAcid;
import bioinfo.proteins.AminoAcidName;
import bioinfo.proteins.PDBEntry;
import bioinfo.proteins.PDBFileReader;

/**
 * shuffles sequences of PDBEntries
 * will be needed to benchmark any rescoring methods
 * @author andreseitz
 *
 */
public class PDBSequenceShuffle {
	
	private final Random rand = new Random(System.currentTimeMillis());
	/**
	 * bgDistribution contains info to classifiy values from rand into right aminoacids
	 */
	private final double[] bgDistributionBorders = new double[26];
	
	/**
	 * sets bgDistribution in a manner that every aminoacid will be equal distributed
	 */
	private void setStdBgDistributionBorders(){
		double[] tmp = new double[26];
		char[] validAA = {'A','R','N','D','C','E','Q','G','H','I','L','K','M','F','P','S','T','W','Y','V'};
		for(char aa : validAA){
			tmp[aa-65] = (1.0d/20.0d);
		}
		bgDistributionBorders[0] = tmp[0];
		for(int i = 1; i != tmp.length; i++){
			bgDistributionBorders[i] = bgDistributionBorders[i-1]+tmp[i];
		}
		if(bgDistributionBorders[bgDistributionBorders.length] != 1.0d){
			bgDistributionBorders[bgDistributionBorders.length] = 1.0d;
		}
	}
	
	/**
	 * 
	 * @return next random AminoAcidName dependent on given distribution and java.util.Random initialized with current time as seed
	 */
	private AminoAcidName getNextRandomAminoAcidName(){
		double r = rand.nextDouble();
		int tmp = 0;
		while(r > bgDistributionBorders[tmp]){
			tmp++;
		}
		return AminoAcidName.getAAFromOLC((char)(tmp+65));
	}
	
	/**
	 * calculates bgDistribution on one model (preferably the same one that will be shuffled later)
	 * @param model
	 */
	private void setBgDistributionBorders(PDBEntry model){
		if(model == null){
			return;
		}
		int[] tmp = new int[26];
		int index;
		int sum = 0;
		for(int i = 0; i != model.length(); i++){
			index = model.getAminoAcid(i).getName().getOneLetterCode().charAt(0)-65;
			if(index != 'U'-65){
				tmp[index]++;
				sum++;
			}
		}
		bgDistributionBorders[0] = (1.0d*tmp[0]/sum);
		for(int i = 1; i != tmp.length; i++){
			bgDistributionBorders[i] = bgDistributionBorders[i-1]+(1.0d*tmp[i]/sum);
		}
	}
	
	/**
	 * calculates bgDistribution on a set of models
	 * @param models
	 */
	private void setBgDistributionBorders(List<PDBEntry> models){
		int[] tmp = new int[26];
		int index;
		int sum = 0;
		PDBEntry model;
		
		for(int m = 0; m != models.size(); m++){
			model = models.get(m);
			for(int i = 0; i != model.length(); i++){
				index = model.getAminoAcid(i).getName().getOneLetterCode().charAt(0)-65;
				if(index != 'U'-65){
					tmp[index]++;
					sum++;
				}
			}
		}
		bgDistributionBorders[0] = (1.0d*tmp[0]/sum);
		for(int i = 1; i != tmp.length; i++){
			bgDistributionBorders[i] = bgDistributionBorders[i-1]+(1.0d*tmp[i]/sum);
		}
	}
	
	/**
	 * calculates bgDistribution on a set of models
	 * @param models
	 */
	private void setBgDistributionBorders(PDBEntry[] models){
		int[] tmp = new int[26];
		int index;
		int sum = 0;
		PDBEntry model;
		
		for(int m = 0; m != models.length; m++){
			model = models[m];
			for(int i = 0; i != model.length(); i++){
				index = model.getAminoAcid(i).getName().getOneLetterCode().charAt(0)-65;
				if(index != 'U'-65){
					tmp[index]++;
					sum++;
				}
			}
		}
		bgDistributionBorders[0] = (1.0d*tmp[0]/sum);
		for(int i = 1; i != tmp.length; i++){
			bgDistributionBorders[i] = bgDistributionBorders[i-1]+(1.0d*tmp[i]/sum);
		}
	}
	
	/**
	 * calculates bgDistribution on a set of models by idlist and a pdbfilreader with initialized folder
	 * @param models
	 */
	private void setBgDistributionBorders(PDBFileReader reader, List<String> models){
		int[] tmp = new int[26];
		int index;
		int sum = 0;
		PDBEntry model;
		
		for(int m = 0; m != models.size(); m++){
			model = reader.readFromFolderById(models.get(m));
			for(int i = 0; i != model.length(); i++){
				index = model.getAminoAcid(i).getName().getOneLetterCode().charAt(0)-65;
				if(index != 'U'-65){
					tmp[index]++;
					sum++;
				}
			}
		}
		bgDistributionBorders[0] = (1.0d*tmp[0]/sum);
		for(int i = 1; i != tmp.length; i++){
			bgDistributionBorders[i] = bgDistributionBorders[i-1]+(1.0d*tmp[i]/sum);
		}
	}
	
	/**
	 * initializes Object with equal distributed AminoAcidNames
	 * @param models
	 */
	public PDBSequenceShuffle(){
		this.setStdBgDistributionBorders();
	}
	
	/**
	 * initializes Object with equal distributed AminoAcidNames or empty Distribution
	 * CAVE with param false and no further set of the bgDistribution, every shuffle will fail
	 * @param models
	 */
	public PDBSequenceShuffle(boolean setStandard){
		if(setStandard){
			this.setStdBgDistributionBorders();
		}
	}
	
	/**
	 * initializes object with bgDistribution derived from one model
	 * @param model
	 */
	public PDBSequenceShuffle(PDBEntry model){
		this.setBgDistributionBorders(model);
	}
	
	/**
	 * initializes object with bgDistribution derived from a set of models
	 * @param model
	 */
	public PDBSequenceShuffle(List<PDBEntry> models){
		this.setBgDistributionBorders(models);
	}
	
	/**
	 * initializes object with bgDistribution derived from a set of models
	 * @param model
	 */
	public PDBSequenceShuffle(PDBEntry[] models){
		this.setBgDistributionBorders(models);
	}
	
	/**
	 * initializes object with bgDistribution derived from a set of models
	 * @param model
	 */
	public PDBSequenceShuffle(PDBFileReader reader, List<String> models){
		this.setBgDistributionBorders(reader, models);
	}
	
	/**
	 * shuffles given model with existing distribution
	 * @param model
	 */
	public PDBEntry getShuffled(PDBEntry model){
		return this.getShuffled(model,false);
	}
	
	/**
	 * shuffles given model
	 * @param model
	 * @param updateDistribution will first calculate new distribution if set to true
	 */
	public PDBEntry getShuffled(PDBEntry model, boolean updateDistribution){
		if(updateDistribution){
			this.setBgDistributionBorders(model);
		}
		AminoAcid[] aas = new AminoAcid[model.length()];
		for(int i = 0; i != model.length(); i++){
			aas[i] = new AminoAcid(getNextRandomAminoAcidName(),i, model.getAminoAcid(i).getAtoms());
		}
		return new PDBEntry(model.getID()+model.getChainID()+String.format("%02d",model.getChainIDNum())+"_shuffled", aas);
	}
	
	/**
	 * shuffles given models with existing distribution
	 * @param models
	 */
	public PDBEntry[] getShuffled(PDBEntry[] models){
		return this.getShuffled(models,false);
	}
	
	/**
	 * shuffles given models
	 * @param model
	 * @param updateDistribution will first calculate new distribution if set to true
	 */
	public PDBEntry[] getShuffled(PDBEntry[] models, boolean updateDistribution){
		if(updateDistribution){
			this.setBgDistributionBorders(models);
		}
		AminoAcid[] aas;
		PDBEntry[] results = new PDBEntry[models.length];
		PDBEntry model;
		for(int m = 0; m != models.length; m++){
			model = models[m];
			aas = new AminoAcid[model.length()];
			for(int i = 0; i != model.length(); i++){
				aas[i] = new AminoAcid(getNextRandomAminoAcidName(),i, model.getAminoAcid(i).getAtoms());
			}
			results[m] = new PDBEntry(model.getID()+model.getChainID()+String.format("%02d",model.getChainIDNum())+"_shuffled"+m, aas);
		}
		return results;
	}
	
	/**
	 * shuffles given models with existing distribution
	 * @param models
	 */
	public List<PDBEntry> getShuffled(List<PDBEntry> models){
		return this.getShuffled(models,false);
	}
	
	/**
	 * shuffles given models
	 * @param model
	 * @param updateDistribution will first calculate new distribution if set to true
	 */
	public List<PDBEntry> getShuffled(List<PDBEntry> models, boolean updateDistribution){
		if(updateDistribution){
			this.setBgDistributionBorders(models);
		}
		AminoAcid[] aas;
		List<PDBEntry> results = new ArrayList<PDBEntry>();
		PDBEntry model;
		for(int m = 0; m != models.size(); m++){
			model = models.get(m);
			aas = new AminoAcid[model.length()];
			for(int i = 0; i != model.length(); i++){
				aas[i] = new AminoAcid(getNextRandomAminoAcidName(),i, model.getAminoAcid(i).getAtoms());
			}
			results.add(new PDBEntry(model.getID()+model.getChainID()+String.format("%02d",model.getChainIDNum())+"_shuffled"+m, aas));
		}
		return results;
	}
	
	/**
	 * generates random sequence on given length and distribution previously set to shuffle instance.
	 * @param length
	 * @return random amino acid name sequence
	 */
	public AminoAcidName[] shuffleSequence(int length){
		AminoAcidName[] result = new AminoAcidName[length];
		for(int i = 0; i != length; i++){
			result[i] = getNextRandomAminoAcidName();
		}
		return result;
	}

}
