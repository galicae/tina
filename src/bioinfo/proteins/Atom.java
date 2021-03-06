/******************************************************************************
 * bioinfo.proteins.Atom                                                      *
 *                                                                            *
 * This file is best read at line width 80 and tab width 4.                   *
 ******************************************************************************/
package bioinfo.proteins;

import java.io.Serializable;
import java.text.DecimalFormat;

/**
 * @author gobi_4
 * @date November 23, 2012
 */
public class Atom implements Serializable{

	/**
	 * 
	 */
	private static final long serialVersionUID = -3326166733645976909L;
	
	private final AtomType type;
	private double[] position;

	public Atom(AtomType type, double[] position) {
		this.type = type;
		this.position = position;
	}

	public Atom(String type, double[] position) {
		this.type = AtomType.createFromString(type);
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
		return type + "\t" + position[0] + "\t" + position[1] + "\t"
				+ position[2] + "\n";
	}

	/**
	 * this function returns a line representation of the atom
	 * 
	 * @param atomNumber
	 *            the number of the atom
	 * @param resNumber
	 *            the number of the amino acid to which the atom belongs
	 * @param resName
	 *            the name of the amino acid to which the atom belongs
	 * @param chain
	 *            the chain to which the atom belongs
	 * @return a string representation of the string that should be parsable by
	 *         any PDB reader
	 */
	public String toString(int atomNumber, int resNumber, String resName,
			char chain) {
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
		// 1 2 3 4 5 6 7 8
		// 12345678901234567890123456789012345678901234567890123456789012345678901234567890
		// ATOM 296 N ALA A 1 10.883 6.779 -6.464 1.00 0.00 N
		String result = "ATOM  ##### aaaa+rrr c****i   xxxxxxxxyyyyyyyyzzzzzzzzooooootttttt          eehh";

		// coordinate strings
		DecimalFormat df = new DecimalFormat("###0.000");
		String xCoord = df.format(position[0]).replace(",", ".");
		while (xCoord.length() < 8)
			xCoord = " " + xCoord;
		String yCoord = df.format(position[1]).replace(",", ".");
		while (yCoord.length() < 8)
			yCoord = " " + yCoord;
		String zCoord = df.format(position[2]).replace(",", ".");
		while (zCoord.length() < 8)
			zCoord = " " + zCoord;
		// atom type
		String atomType = type.toString();
		if (atomType.length() < 4)
			atomType = " " + atomType;
		while (atomType.length() < 4)
			atomType = atomType + " ";
		// element symbol
		String element = "";
		if (atomType.contains("C"))
			element = "C";
		else if (atomType.contains("N"))
			element = "N";
		else if (atomType.contains("O"))
			element = "O";
		else
			element = Character.toString(atomType.charAt(0));
		while (element.length() < 2)
			element = " " + element;
		
		if(resName.length() == 1) {
			resName = AminoAcidName.getThreeLetterCode(resName);
		}
		result = result.replace("oooooo", "  1.00");
		result = result.replace("tttttt", "  0.00");
		result = result.replace("h", " ");
		result = result.replace("c", chain + "");
		result = result.replace("+", " ");
		result = result.replace("i", " ");
		result = result.replace("ee", element);
		result = result.replace("****", String.format("%4d", resNumber));
		result = result.replace("rrr", resName);
		result = result.replace("#####", String.format("%5d", atomNumber));
		result = result.replace("aaaa", atomType);
		result = result.replace("xxxxxxxx", xCoord);
		result = result.replace("yyyyyyyy", yCoord);
		result = result.replace("zzzzzzzz", zCoord);
		return result;
	}

	
	/**
	 * this function is the same as the toString function, although without the
	 * extra information. It uses placeholders instead, making the creation of a
	 * line of length 80 where everything is in its place easier
	 * 
	 * @return
	 */
	public String toTMString() {
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
		// 1 2 3 4 5 6 7 8
		// 12345678901234567890123456789012345678901234567890123456789012345678901234567890
		// ATOM 296 N ALA A 1 10.883 6.779 -6.464 1.00 0.00 N
		String result = "ATOM  ##### aaaa+rrr c****i   xxxxxxxxyyyyyyyyzzzzzzzzooooootttttt          eehh";

		// coordinate strings
		String xCoord = Double.toString(position[0]);
		while (xCoord.length() < 8)
			xCoord = " " + xCoord;
		String yCoord = Double.toString(position[1]);
		while (yCoord.length() < 8)
			yCoord = " " + yCoord;
		String zCoord = Double.toString(position[2]);
		while (zCoord.length() < 8)
			zCoord = " " + zCoord;
		// atom type
		String atomType = type.toString();
		while (atomType.length() < 4)
			atomType = " " + atomType;
		// element symbol
		String element = "";
		if (atomType.contains("C"))
			element = "C";
		else if (atomType.contains("N"))
			element = "N";
		else if (atomType.contains("O"))
			element = "O";
		else
			element = Character.toString(atomType.charAt(0));
		while (element.length() < 2)
			element = " " + element;

		// result = result.replace("o", " ");
		// result = result.replace("t", " ");
		// result = result.replace("h", " ");
		// result = result.replace("c", chain+"");
		// result = result.replace("+", " ");
		// result = result.replace("i", " ");
		result = result.replace("ee", element);
		// result = result.replace("****", String.format("%4d",resNumber));
		// result = result.replace("rrr", resName);
		// result = result.replace("#####", String.format("%5d",atomNumber));
		result = result.replace("aaaa", atomType);
		result = result.replace("xxxxxxxx", xCoord);
		result = result.replace("yyyyyyyy", yCoord);
		result = result.replace("zzzzzzzz", zCoord);
		return result;
	}
	
	public Atom clone() {
		Atom clone = new Atom(this.type, this.position);
		return clone;
	}
}
