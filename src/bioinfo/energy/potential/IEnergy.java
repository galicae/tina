package bioinfo.energy.potential;

import java.util.List;
import java.util.Map;

public interface IEnergy {
	
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
	 * calculates energy from PDB in a Boltzmann-based mean force way
	 * should only have to be done once, written to file and then reused from there
	 * @param List of pdb ids, which should be used to generate energy function
	 * location of PDBs can be implemented in abstract class
	 */
	public void calculateFromPDBs(List<String> pdbIds);
	
	/**
	 * 
	 * @param mapping parameter set, passed as key - value pairs
	 * @return double value containing energy value for given set of Paramters in mapping
	 */
	public double getEnergyValue(Map<Object, Object> mapping);
	
	/**
	 * Datatypes of values to given keys should be taken from the method in Energy Potential implementing this interface
	 * @return List of String values representing valid (and required!) keys which can be used with getEnergyValue
	 */
	public String[] getMapKeys();
	
}