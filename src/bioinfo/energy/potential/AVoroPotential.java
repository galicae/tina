package bioinfo.energy.potential;

import bioinfo.energy.potential.voronoi.VoroPPWrap;
import bioinfo.energy.potential.voronoi.VoroPrepType;
import bioinfo.energy.potential.voronoi.VoronoiData;
import bioinfo.proteins.DSSPEntry;
import bioinfo.proteins.PDBEntry;

public abstract class AVoroPotential extends APotential{
	
	protected VoronoiData data;
	protected VoroPPWrap voro;
	
	public AVoroPotential(String vorobin, String tmpdir){
		voro = new VoroPPWrap(tmpdir,vorobin);
	}
	
	public AVoroPotential(String vorobin){
		voro = new VoroPPWrap("/tmp/", vorobin);
	}
	
	public VoronoiData prepareSimple(PDBEntry pdb, VoroPrepType type){
		data = new VoronoiData(pdb.getID());
		data.reducePDB(type, pdb);
		voro.decomposite(data);
		return data;
	}
	
	public VoronoiData prepareSimple(DSSPEntry dssp){
		data = new VoronoiData(dssp.getID());
		data.reduceDSSP(dssp);
		voro.decomposite(data);
		return data;
	}
	
	public VoronoiData prepareWithGrid(PDBEntry pdb, double gridExtend, double gridDensity, double gridClash, double minContact){
		data = new VoronoiData(pdb.getID());
		data.reducePDB(pdb);
		data.fillGridWithoutClashes(gridExtend, gridDensity, gridClash);
		voro.decomposite(data);
		data.detectOuterGrid(minContact);
		return data;
	}
	
	public VoronoiData prepareWithGrid(DSSPEntry dssp, double gridExtend, double gridDensity, double gridClash, double minContact){
		data = new VoronoiData(dssp.getID());
		data.reduceDSSP(dssp);
		data.fillGridWithoutClashes(gridExtend, gridDensity, gridClash);
		voro.decomposite(data);
		data.detectOuterGrid(minContact);
		return data;
	}		

}
