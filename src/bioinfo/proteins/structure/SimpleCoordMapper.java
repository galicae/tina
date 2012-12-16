/**
 * 
 */
package bioinfo.proteins.structure;

import java.util.LinkedList;

import bioinfo.alignment.SequenceAlignment;
import bioinfo.proteins.AminoAcid;
import bioinfo.proteins.Atom;
import bioinfo.proteins.AtomType;
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
				//AminoAcid temp = arg2.getAminoAcid(map[0][i]);
				// TODO debugging: is the PDBEntry given at all?
				if (arg2 == null) {
					System.err.println("debugging: The given PDBEntry is null!");
				}
				if (arg2.getAminoAcid(map[0][i]) != null) {
					// DONE: only map correct CA or Backbone or fitting AAs
					// momentarily it only maps the backbone
					AminoAcid fromtemp = arg2.getAminoAcid(map[0][i]);
					LinkedList<Atom> tempAtoms = new LinkedList<Atom>();
					
					// N
					Atom tempatom = fromtemp.getAtomByType(AtomType.N);
					if (tempatom != null) {
						tempAtoms.add(tempatom);
					}
					// CA
					tempatom = fromtemp.getAtomByType(AtomType.CA);
					if (tempatom != null) {
						tempAtoms.add(tempatom);
					}
					// C
					tempatom = fromtemp.getAtomByType(AtomType.C);
					if (tempatom != null) {
						tempAtoms.add(tempatom);
					}
					// O
					tempatom = fromtemp.getAtomByType(AtomType.O);
					if (tempatom != null) {
						tempAtoms.add(tempatom);
					}
					
					AminoAcid temp = new AminoAcid(arg2.getAminoAcid(map[0][i]).getName(), tempAtoms.toArray(new Atom[0]));
					aalist.add(temp);
				}
			}
		}
		
		return new PDBEntry(arg1.getComponent(0).getID(), aalist);
	}
	
}
