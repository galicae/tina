package bioinfo.energy.potential.preparation.voronoi;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

import cern.colt.matrix.DoubleFactory1D;

import bioinfo.proteins.AminoAcid;
import bioinfo.proteins.Atom;
import bioinfo.proteins.AtomType;
import bioinfo.proteins.PDBEntry;
import bioinfo.proteins.PDBFileReader;

import opt.CommandLine;
import opt.ParameterType;
import opt.Setting;
import opt.ValidationRule;


/**
 * 
 * @author andreseitz
 *	VoroPrep prepares PDB-Files in the following manner:
 *	output will be a \t separated file containing <numerical-id>,xcoord,ycoord,zcoord
 *	different paramters have to be passed like type of coordinate and pdbfile
 *
 */
public class VoroPrepare {
	
	/**
	 * @param args
	 */
	public VoronoiData reducePDB(VoroPrepType type, PDBEntry pdb){

		VoronoiData reduced = new VoronoiData(pdb.getID(),type);
		
		switch(type){
			case CA:{
				Atom tmp = null;
				for(int i = 0; i != pdb.length(); i++){
					if((tmp = pdb.getAminoAcid(i).getAtomByType(AtomType.CA)) != null){
						reduced.addValues(i,tmp.getPosition(),pdb.getAminoAcid(i).getName());
					}
				}
			};break;
			/*
			case AA:{
				AminoAcid tmp;
				Atom temp;
				int counter = 0;
				for(int i = 0; i != pdb.length(); i++){
					tmp = pdb.getAminoAcid(i);
					for(int j = 0; j != tmp.getAtomNumber(); j++){
						temp = tmp.getAtom(j);
						reduced.addValues(counter, temp.getPosition(), tmp.getName())
						counter++;
					}
				}
			};break;
			*/
			case CC:{
				AminoAcid tmp;
				Atom temp;
				double[][] coord;
				double[] centroid;
				for(int i = 0; i != pdb.length(); i++){
					tmp = pdb.getAminoAcid(i);
					coord = new double[tmp.getAtomNumber()][3];
					for(int j = 0; j != tmp.getAtomNumber(); j++){
						temp = tmp.getAtom(j);
						coord[j] = temp.getPosition();
					}
					centroid = calculateCentroid(coord);
					reduced.addValues(i, centroid, tmp.getName());				
				}
			};break;
		}
		
		return reduced;

	}
	
	private static double[] calculateCentroid(double[][] coord){
		double cx = 0, cy = 0, cz = 0;
		for (int i = 0; i < coord.length; i++) {
			cx += coord[i][0];
			cy += coord[i][1];
			cz += coord[i][2];
		}
		cx = cx / (coord.length * 1.0);
		cy = cy / (coord.length * 1.0);
		cz = cz / (coord.length * 1.0);
		double[] result = { cx, cy, cz };
		return result;
	}

}
