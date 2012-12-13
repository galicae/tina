/**
 * 
 */
package bioinfo.proteins.structure;

import java.util.LinkedList;

import bioinfo.alignment.SequenceAlignment;
import bioinfo.proteins.AminoAcid;
import bioinfo.proteins.PDBEntry;

/**
 * @author huberste
 *
 */
public class SimpleCoordMapper {
	
	public static PDBEntry map(SequenceAlignment arg1, PDBEntry arg2) {
		
		int[][] map = arg1.calcMap();
		LinkedList<AminoAcid> aalist = new LinkedList<AminoAcid>();
		
		for (int i = 0; i < map[0].length; i++) {
			if (map[0][i] != -1) { // this aminoacid is aligned!
				aalist.add(arg2.getAminoAcid(map[0][i]).copy());
			}
		}
		
		return new PDBEntry(arg1.getComponent(0).getID(), aalist);
	}
	
}
