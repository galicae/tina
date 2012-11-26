/******************************************************************************
 * bioinfo.proteins.structure.PDBMapper                                       *
 *                                                                            *
 * This file is best read at line width 80 and tab width 4.                   *
 ******************************************************************************/
package bioinfo.proteins.structure;

import bioinfo.alignment.Alignment;
import bioinfo.proteins.PDBEntry;

/**
 * @author gobi_4
 * @date November 24, 2012
 */
public abstract class PDBMapper {

	public PDBMapper() {
		
	}
	
	public abstract PDBEntry map (Alignment arg1, PDBEntry arg2);
	
}