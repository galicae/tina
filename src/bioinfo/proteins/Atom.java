/******************************************************************************
 * bioinfo.proteins.Atom                                                      *
 *                                                                            *
 * This file is best read at line width 80 and tab width 4.                   *
 ******************************************************************************/
package bioinfo.proteins;

import java.util.zip.DataFormatException;

/**
 * @author gobi_4
 * @date November 23, 2012
 */
public class Atom {

	private final AtomType type;
	private double[] position;

	public Atom(AtomType type, double[] position) {
		this.type = type;
		this.position = position;
	}

	public Atom(String type, double[] position) {
		if (type.trim().equals("CA"))
			this.type = AtomType.CA;
		else if (type.trim().equals("C"))
			this.type = AtomType.C;
		else if (type.trim().equals("N"))
			this.type = AtomType.N;
		else if (type.trim().equals("O"))
			this.type = AtomType.O;
		else if (type.trim().equals("CB"))
			this.type = AtomType.CB;
		else if (type.trim().equals("CG"))
			this.type = AtomType.CG;
		else if (type.trim().equals("CD"))
			this.type = AtomType.CD;
		else if (type.trim().equals("NE"))
			this.type = AtomType.NE;
		else if (type.trim().equals("CZ"))
			this.type = AtomType.CZ;
		else if (type.trim().equals("NH1"))
			this.type = AtomType.NH1;
		else if (type.trim().equals("NH2"))
			this.type = AtomType.NH2;
		else
			this.type = AtomType.U;
		this.position = position;
	}

	public AtomType getType() {
		return type;
	}

	public double[] getPosition() {
		return position;
	}

	public void setPosition(double[] array) {
		position = array;
	}

	@Override
	public String toString() {
		// how to produce perfect PDB atom records:
		// #: atom serial number
		// a: atom name
		// +: alternate location indicator
		// r: residue name
		// c: chain identifier
		// *: residue sequence number
		// i: code for insertion of residues
		// x: x coordinate (8.3) // 8 length overall, 3 digits after point
		// y: y coordinate (8.3)
		// z: z coordinate (8.3)
		// o: occupancy (6.2)
		// t: temperature factor (6.2)
		// e: element symbol
		// h: charge of the atom
		// element symbol (e) should always be right-justified
//                                1         2         3         4         5         6         7         8
//                       12345678901234567890123456789012345678901234567890123456789012345678901234567890
        //               ATOM    296  N   ALA  A   1     10.883   6.779  -6.464  1.00  0.00           N
		String result = "ATOM  ##### aaaa+rrr c****i   xxxxxxxxyyyyyyyyzzzzzzzzooooootttttt          eehh";
		
//		coordinate strings
		String xCoord = Double.toString(position[0]);
		while(xCoord.length() < 8)
			xCoord = " " + xCoord;
		String yCoord = Double.toString(position[0]);
		while(yCoord.length() < 8)
			yCoord = " " + yCoord;
		String zCoord = Double.toString(position[0]);
		while(zCoord.length() < 8)
			zCoord = " " + zCoord;
//		atom type
		String atomType = type.toString();
		while(atomType.length() < 4)
			atomType = " " + atomType;
//		element symbol
		String element = "";
		if(atomType.contains("C"))
			element = "C";
		else if(atomType.contains("N"))
			element = "N";
		else if(atomType.contains("O"))
			element = "O";
		else
			element = Character.toString(atomType.charAt(0));
		while(element.length() < 2)
			element = " " + element;
		
		result.replace("ee", element);
		result.replace("aaaa", atomType);
		result.replace("xxxxxxxx", xCoord);
		result.replace("yyyyyyyy", yCoord);
		result.replace("zzzzzzzz", zCoord);
		return result;
	}
}
