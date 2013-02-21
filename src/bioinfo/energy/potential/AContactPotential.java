package bioinfo.energy.potential;

import bioinfo.proteins.PDBEntry;

public abstract class AContactPotential extends APotential{
	
	/**
	 * model only needed and set for native vs shuffle scoring with potential
	 */
	protected PDBEntry model;
	
	/**
	 * initializes bulk readout of shuffled sequences on model
	 * @param model
	 */
	public void prepareSequenceScoring(Object model){
		this.model = (PDBEntry)model;
	}
}
