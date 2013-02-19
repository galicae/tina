package bioinfo.energy.potential;

import bioinfo.energy.potential.voronoi.VoronoiData;
import bioinfo.proteins.PDBEntry;

public abstract class AContactPotential extends APotential{
	
	/**
	 * model only needed and set for native vs shuffle scoring with potential
	 */
	protected PDBEntry model;
	
	public void prepareSequenceScoring(PDBEntry model){
		this.model = model;
	}
}
