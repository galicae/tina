package bioinfo.superpos;

import bioinfo.proteins.PDBEntry;
import cern.colt.matrix.DoubleMatrix1D;
import cern.colt.matrix.DoubleMatrix2D;

/**
 * value-object class for transformations which are the output of kabsch
 * @author gobi4
 *
 */

public class Transformation {
	
	private DoubleMatrix1D centroidP;
	private DoubleMatrix1D centroidQ;
	private DoubleMatrix2D rotation;
	private double rmsd;
	
	/**
	 * 
	 * @param Q PDBEntry which was used for kabsch calculation an has to be transformed on P now
	 * @return new PDBEntry Q which is transformed old Q
	 */
	public PDBEntry transform (PDBEntry Q){
		//TODO
		return null;
	}

}
