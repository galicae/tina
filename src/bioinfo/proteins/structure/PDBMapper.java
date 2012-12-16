/******************************************************************************
 * bioinfo.proteins.structure.PDBMapper                                       *
 *                                                                            *
 * This file is best read at line width 80 and tab width 4.                   *
 ******************************************************************************/
package bioinfo.proteins.structure;

import bioinfo.alignment.SequenceAlignment;
import bioinfo.proteins.PDBEntry;

/**
 * @author gobi_4
 * @date November 24, 2012
 */
<<<<<<< HEAD
public abstract class PDBMapper {

	public PDBMapper(){
		
	}
=======
public interface PDBMapper {
>>>>>>> branch 'master' of https://github.com/galicae/tina.git
	
	public PDBEntry map (SequenceAlignment arg1, PDBEntry arg2);
	
}
