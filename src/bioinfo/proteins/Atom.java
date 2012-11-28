/******************************************************************************
 * bioinfo.proteins.Atom                                                      *
 *                                                                            *
 * This file is best read at line width 80 and tab width 4.                   *
 ******************************************************************************/
package bioinfo.proteins;

/**
 * @author gobi_4
 * @date November 23, 2012
 */
public class Atom {

	
	private final AtomType type;
	private final double[] position;
	
	public Atom(AtomType type, double[] position) {
		this.type = type;
		this.position = position;
	}
	
	public Atom(String type, double[] position) {
		if (type.trim().equals("CA")) this.type = AtomType.CA;
		else if (type.trim().equals("C")) this.type = AtomType.C;
		else if (type.trim().equals("N")) this.type = AtomType.N;
		else if (type.trim().equals("O")) this.type = AtomType.O;
		else if (type.trim().equals("CB")) this.type = AtomType.CB;
		else if (type.trim().equals("CG")) this.type = AtomType.CG;
		else if (type.trim().equals("CD")) this.type = AtomType.CD;
		else if (type.trim().equals("NE")) this.type = AtomType.NE;
		else if (type.trim().equals("CZ")) this.type = AtomType.CZ;
		else if (type.trim().equals("NH1")) this.type = AtomType.NH1;
		else if (type.trim().equals("NH2")) this.type = AtomType.NH2;
		else this.type = AtomType.U;
		this.position = position;
	}
	
	public AtomType getType() {
		return type;
	}
	
	public double[] getPosition() {
		return position;
	}
	
}
