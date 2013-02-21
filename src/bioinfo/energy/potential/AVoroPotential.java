package bioinfo.energy.potential;

import java.util.Set;

import bioinfo.energy.potential.voronoi.VoroPPWrap;
import bioinfo.energy.potential.voronoi.VoroPrepType;
import bioinfo.energy.potential.voronoi.VoronoiData;
import bioinfo.proteins.AminoAcidName;
import bioinfo.proteins.DSSPEntry;
import bioinfo.proteins.PDBEntry;

public abstract class AVoroPotential extends APotential{
	
	/**
	 * holds voronoi wrapper
	 *
	 */
	protected VoroPPWrap voro;
	
	/**
	 * data container for voronoi decomposited models
	 * only generated when shuffledSequence scoring shall be used
	 */
	protected VoronoiData data;
	
	/**
	 * data containers for voronoi decomposited parameters
	 * only generated when shuffledSequence scoring shall be used
	 */
	protected double minContact;
	protected double gridExtend;
	protected double gridDensity;
	protected double gridClash;
	
	/**
	 * @param vorobin binary position
	 * @param tmpdir eg /tmp/
	 */
	public AVoroPotential(String vorobin, String tmpdir, double minContact, double gridExtend, double gridDensity, double gridClash){
		voro = new VoroPPWrap(tmpdir,vorobin);
		this.minContact = minContact;
		this.gridExtend = gridExtend;
		this.gridDensity = gridDensity;
		this.gridClash = gridClash;
	}
	
	/**
	 * @param voro params as 4 double values
	 * should only be used after initializing shuffle sequence read out of the same structure
	 */
	public void setVoroParam(String vorobin, double minContact, double gridExtend, double gridDensity, double gridClash){
		this.minContact = minContact;
		this.gridExtend = gridExtend;
		this.gridDensity = gridDensity;
		this.gridClash = gridClash;
	}
	
	/**
	 * @param vorobin binary position
	 * @param tmpdir eg /tmp/
	 */
	public AVoroPotential(String vorobin, String tmpdir){
		voro = new VoroPPWrap(tmpdir,vorobin);
	}
	
	/**
	 * @param vorobin binary position
	 * tempdir will be set to /tmp/
	 */
	public AVoroPotential(String vorobin){
		voro = new VoroPPWrap("/tmp/", vorobin);
	}
	
	/**
	 * 
	 * @param pdb
	 * @param type
	 * @return voronoi data reuced and decomposited by given pdbmodel
	 * no grid will be used
	 */
	public VoronoiData prepareSimple(PDBEntry pdb, VoroPrepType type){
		VoronoiData data;
		data = new VoronoiData(pdb.getId());
		data.reducePDB(type, pdb);
		voro.decomposite(data);
		return data;
	}
	
	/**
	 * 
	 * @param dssp
	 * @return voronoi data reuced and decomposited by given pdbmodel
	 * no grid will be used!
	 */
	public VoronoiData prepareSimple(DSSPEntry dssp){
		VoronoiData data;
		data = new VoronoiData(dssp.getId());
		data.reduceDSSP(dssp);
		voro.decomposite(data);
		return data;
	}
	
	/**
	 * normal preparation of data in solvent
	 * @param pdb
	 * @param gridExtend
	 * @param gridDensity
	 * @param gridClash
	 * @param minContact
	 * @return voronoi data reuced and decomposited by given pdbmodel
	 * grid will be used and generated
	 */
	public VoronoiData prepareWithGrid(PDBEntry pdb, double gridExtend, double gridDensity, double gridClash, double minContact){
		VoronoiData data;
		data = new VoronoiData(pdb.getId());
		data.reducePDB(pdb);
		data.fillGridWithoutClashes(gridExtend, gridDensity, gridClash);
		voro.decomposite(data);
		data.detectOuterGrid(minContact);
		return data;
	}
	
	/**
	 * normal preparation of data in solvent
	 * @param dssp
	 * @param gridExtend
	 * @param gridDensity
	 * @param gridClash
	 * @param minContact
	 * @return voronoi data reuced and decomposited by given pdbmodel
	 * grid will be generated and used
	 */
	public VoronoiData prepareWithGrid(DSSPEntry dssp, double gridExtend, double gridDensity, double gridClash, double minContact){
		VoronoiData data;
		data = new VoronoiData(dssp.getId());
		data.reduceDSSP(dssp);
		data.fillGridWithoutClashes(gridExtend, gridDensity, gridClash);
		voro.decomposite(data);
		data.detectOuterGrid(minContact);
		return data;
	}		
	
	/**
	 * voronoi data container for multiple sequence scoring on model will be initialized
	 * @param model
	 * @return voronoi data reuced and decomposited by given pdbmodel
	 * grid will be used and generated
	 */
	@Override
	public void prepareSequenceScoring(Object model){
		PDBEntry entry = (PDBEntry)model;
		data = new VoronoiData(entry.getId());
		data.reducePDB(entry);
		data.fillGridWithoutClashes(gridExtend, gridDensity, gridClash);
		voro.decomposite(data);
		data.detectOuterGrid(minContact);
	}
	
	/**
	 * return s amino acid wise scores in double array
	 * no decomposition shall be done in this method
	 * @param data
	 * @return 
	 */
	protected abstract double[] getAminoScores(VoronoiData data);
	
	/**
	 * returns score of native model set previously by prepareSquenceScoring
	 * no more decomposition will be made by calling this method
	 */
	@Override
	public double getNativeScoring(){
		double[] scores = getNativeAminoScoring();
		double score = 0.0d;
		for(int i = 0; i != scores.length; i++){
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

	
	

}
