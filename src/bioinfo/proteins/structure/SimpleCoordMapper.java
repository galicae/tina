/******************************************************************************
 * bioinfo.proteins.structure.SimpleCoordMapper                               *
 *                                                                            *
 * This file is best read at line width 80 and tab width 4.                   *
 *                                                                   huberste *
 ******************************************************************************/
package bioinfo.proteins.structure;

import java.util.LinkedList;

import bioinfo.alignment.SequenceAlignment;
import bioinfo.proteins.AminoAcid;
import bioinfo.proteins.Atom;
import bioinfo.proteins.AtomType;
import bioinfo.proteins.PDBEntry;

/**
 * SimpleCoordMapper provides a static function to map an PDBEntry onto an
 * alignment.
 * @author huberste
 * @lastchange 2013-02-13
 */
public class SimpleCoordMapper {

	/**
	 * maps an PDBEntry onto an SequenceAlignment
	 * @param arg1 a SequenceAlignment with row[0] the target and row[1] the template.
	 * @param arg2 the PDBEntry for the template
	 * @return the PDBEntry for the target
	 */
	public static PDBEntry map(SequenceAlignment arg1, PDBEntry arg2) throws NullPointerException {

		// debugging: is the PDBEntry given at all?
		if (arg2 == null) {
//			System.err.println(
//				"debugging: The given PDBEntry is null!"
//			);
			throw new NullPointerException();
		}
		
		int[][] map = arg1.calcMap();
		LinkedList<AminoAcid> aalist = new LinkedList<AminoAcid>();

		for (int i = 0; i < map[0].length; i++) {
			if (map[0][i] != -1) { // this aminoacid is aligned!
				// AminoAcid temp = arg2.getAminoAcid(map[0][i]);
				
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
					
					AminoAcid temp = new AminoAcid(arg2.getAminoAcid(map[0][i])
							.getName(), arg2.getAminoAcid(map[0][i])
							.getResIndex(), tempAtoms.toArray(new Atom[0]));
					aalist.add(temp);
				}
			}
		}

		return new PDBEntry(arg1.getComponent(0).getId(), aalist);
	}
	
	/**
	 * maps an PDBEntry onto an SequenceAlignment
	 * @param arg1 a SequenceAlignment with row[0] the template and row[1] the target.
	 * @param arg2 the PDBEntry for the template
	 * @return the PDBEntry for the target
	 */
	public static PDBEntry map(PDBEntry arg2, SequenceAlignment arg1) throws NullPointerException {

		// debugging: is the PDBEntry given at all?
		if (arg2 == null) {
//			System.err.println(
//				"debugging: The given PDBEntry is null!"
//			);
			throw new NullPointerException();
		}
		
		int[][] map = arg1.calcMap();
		LinkedList<AminoAcid> aalist = new LinkedList<AminoAcid>();
		
		// i is the i-th position in the second sequence
		for (int i = 0; i < map[1].length; i++) {
			aalist.add(new AminoAcid(arg1.getComponent(1).getComp(i),i));
			if (map[1][i] != -1) { // this aminoacid in target is aligned!
				
				if (arg2.getAminoAcid(map[1][i]) != null) {
					// DONE: only map correct CA or Backbone or fitting AAs
					// momentarily it only maps the backbone
					AminoAcid fromtemp = arg2.getAminoAcid(map[1][i]);
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
					
					aalist.get(i).setAtoms(tempAtoms.toArray(new Atom[0]));
				}
			}
		}

		return new PDBEntry(arg1.getComponent(1).getId(), aalist);
	}

}
