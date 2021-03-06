package bioinfo.energy.potential;

import java.util.List;
import bioinfo.proteins.AminoAcidName;
import bioinfo.proteins.PDBEntry;

public interface IEnergy {
	
	/**
	 * 
	 * @param size array of size-values for each potential array field
	 * e.g. {2,2} generates potential[2][2]
	 */
	public void initPotential(int[] size);
	
	/**
	 * writes energy to file
	 * @param filename where energy is to be written into
	 */
	public void writeToFile(String filename);
	
	/**
	 * reads energy from file
	 * @param filename where energy was written into
	 */
	public void readFromFile(String filename);
	
	/**
	 * calculates energy from data in a Boltzmann-based mean force way
	 * should only have to be done once, written to file and then reused from there
	 * @param List of ids, which should be used to generate energy function
	 * location of data can be implemented in abstract class
	 * @param dataFolder containing all data to every listed id to calculate
	 * potential from
	 */
	public void calculateFromDATA(List<String> Ids, String dataFolder);
	
	/**
	 * 
	 * @param model of type PDBEntry
	 * @return score of model calculated with potential
	 */
	public double scoreModel(PDBEntry model);
	
	/**
	 * 
	 * @param model of type PDBEntry
	 * @return aminoacid-wise scores of model read from potential
	 */
	public double[] getAminoScores(PDBEntry model);
	
	/**
	 * needed for efficient calculation of evaluation-data
	 * PDBEntry model has to be set in previous step conatining Structure
	 * then one can access native score and a lot of shuffled sequences on 
	 * the same structure without calculating structure dependent infos each time
	 * 
	 * @return score for native sequence on native structure, prepared earlier
	 */
	public double getNativeScoring();
	
	/**
	 * needed for efficient calculation of evaluation-data
	 * PDBEntry model has to be set in previous step conatining Structure
	 * then one can access native score and a lot of shuffled sequences on 
	 * the same structure without calculating structure dependent infos each time
	 * 
	 * @param sequence eg generated by shuffled, that matches length of model
	 * @return score of given sequence on native structure
	 */
	public double getSequenceScoring(AminoAcidName[] sequence);
	
	/**
	 * needed for efficient calculation of evaluation-data
	 * PDBEntry model has to be set in previous step conatining Structure
	 * then one can access native score and a lot of shuffled sequences on 
	 * the same structure without calculating structure dependent infos each time
	 * 
	 * @return footprint scores for native sequence on native structure, prepared earlier
	 */
	public double[] getNativeAminoScoring();
	
	/**
	 * needed for efficient calculation of evaluation-data
	 * PDBEntry model has to be set in previous step conatining Structure
	 * then one can access native score and a lot of shuffled sequences on 
	 * the same structure without calculating structure dependent infos each time
	 * 
	 * @param sequence eg generated by shuffled, that matches length of model
	 * @return footprint scores of given sequence on native structure
	 */
	public double[] getSequenceAminoScoring(AminoAcidName[] sequence);
	
	/**
	 * sets model as native model to be used with shuffled sequences read out as 
	 * bulk without generating decomposition ech time a method is called
	 * @param model to be set as native
	 */
	public void prepareSequenceScoring(Object data);
	
	
}
