package bioinfo.energy.potential;

import java.util.HashMap;
import java.util.List;

import bioinfo.alignment.matrices.QuasarMatrix;
import bioinfo.energy.potential.voronoi.VoronoiData;
import bioinfo.proteins.AminoAcidName;
import bioinfo.proteins.PDBEntry;
import bioinfo.proteins.SecStructThree;
import bioinfo.proteins.structure.PDBSequenceShuffle;

/**
 * Potential class for testing if vorolign diffscore can be used as method for potential generation
 * @author andreseitz
 *
 */
public class VorolignPotential extends AVoroPotential{

	public VorolignPotential(double[][] scoringMatrix){
		super(null);
		this.scoringMatrix = scoringMatrix;
	}
	
	/**
	 * scoringMatrix used for calculating diffscore between "native" and other sequence regarding vorolign contacts
	 */
	double[][] scoringMatrix;
	
	/**
	 * voronoi data container for shuffle sequence bulk readout
	 */
	private VoronoiData data;
	
	/**
	 * nothing to do here, cause Vorloign Potential is no real potential
	 * its scores a calculated on the fly in a vorologn like manner
	 */
	@Override
	public void calculateFromDATA(List<String> Ids, String dataFolder) {
		return;	
	}
	
	/**
	 * voronoi data container for multiple sequence scoring on model will be initialized
	 * @param model
	 * @return voronoi data reuced and decomposited by given pdbmodel
	 * grid will be used and generated
	 */
	@Override
	public void prepareSequenceScoring(Object model){
		data = (VoronoiData)model;
	}

	
	/**
	 * returns score of native model set previously by prepareSquenceScoring
	 * no more decomposition will be made by calling this method
	 */
	@Override
	public double getNativeScoring(){
		double[] scores = getNativeAminoScoring();
		double score = 0.0d;
		for(int i = 0; i!= scores.length; i++){
			score += scores[i];
		}
		return score;
	}
	
	/**
	 * returns amino acid wise score of native model set previously by prepareSquenceScoring
	 * no more decomposition will be made by calling this method
	 */
	@Override
	public double[] getNativeAminoScoring(){
		return getAminoScores(this.data);
	}
	
	/**
	 * returns score of given sequence on model set previously by prepareSquenceScoring
	 * no more decomposition will be made by calling this method
	 */
	@Override
	public double getSequenceScoring(AminoAcidName[] sequence){
		double[] scores = getSequenceAminoScoring(sequence);
		double score = 0.0d;
		for(int i = 0; i != scores.length; i++){
			score += scores[i];
		}
		return score;
	}
	
	/**
	 * returns amino acid wise score of given sequence on model set previously by prepareSquenceScoring
	 * no more decomposition will be made by calling this method
	 */
	@Override
	public double[] getSequenceAminoScoring(AminoAcidName[] sequence){
		VoronoiData data = this.data.clone();
		data.overrideAminoNames(sequence);
		return getAminoScores(data);
	}

	/**
	 * doesnt make sense here
	 * method is used for special investigation of "vorolign alignments"
	 */
	@Override
	public double scoreModel(PDBEntry model) {
		return 0;
	}

	/**
	 * doesnt make sense here
	 * method is used for special investigation of "vorolign alignments"
	 */
	@Override
	public double[] getAminoScores(PDBEntry model) {
		return null;
	}

	/**
	 * prepareSequenceScoring has to be called before!
	 * calculates "diff-score" := matrix(native,native)-matrix(native,sample)
	 */
	@Override
	protected double[] getAminoScores(VoronoiData data) {
		int[] numContacts = new int[data.getPeptideIds().size()];
		HashMap<Integer,Double> faces;
		for(int i = 0; i != numContacts.length; i++){
			faces = this.data.getFaces().get(i);
			for(int f : faces.keySet()){
				if(this.data.getPeptideIds().contains(f)){
					numContacts[i]++;
				}
			}
		}
		double[] scores = new double[data.getPeptideIds().size()];
		int natRes = 0;
		int datRes = 0;
		for(int i = 0; i != scores.length; i++){
			natRes = this.data.getAminos().get(i).getOneLetterCode().charAt(0)-65;
			datRes = data.getAminos().get(i).getOneLetterCode().charAt(0)-65;
			scores[i] = (scoringMatrix[natRes][datRes]-scoringMatrix[natRes][natRes])*numContacts[i];
		}
		return scores;
	}	
	

}
