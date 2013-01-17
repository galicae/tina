/******************************************************************************
 * oracle is the interface for oracles for a rdp.                             *
 *                                                                            *
 * This file is best read at line width 80 and tab width 4.                   *
 *                                                                   huberste *
 ******************************************************************************/
package huberdp;

import bioinfo.proteins.PDBEntry;

/**
 * @author huberste
 * @lastchange 2013-01-17
 */
public interface Oracle {
	
	/**
	 * 
	 * @param template the template structure
	 * @param temp_left left end of template
	 * @param temp_right right end of template
	 * @param target the target structure
	 * @param targ_left left end of target
	 * @param targ_right right end of target
	 * @param number of similiar segments that shall be returned
	 * @return target_structure(s) with only the aligned amino acids
	 */
	public PDBEntry[] findSimiliarSegments(
			PDBEntry template, int temp_left, int temp_right,
			PDBEntry target,   int targ_left, int targ_right,
			int count);
	
}

/******************************************************************************
 * "A question that sometimes drives me hazy:                                 *
 *  Am I or are the others crazy?"                                            *
 *     - Albert Einstein (1979 - 1966)                                        *
 ******************************************************************************/