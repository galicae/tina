/******************************************************************************
 * bioinfo.sequence                                                           *
 *                                                                            *
 * This file is best read at line width 80 and tab width 4.                   *
 ******************************************************************************/
package bioinfo;

import bioinfo.alignment.Alignable;

/**
 * @author gobi_4
 * @date November 25, 2012
 */
public class Sequence implements Alignable{

	private final String id;
	private final char[] sequence;
	private final int length;
	
	public Sequence(String id, char[] sequence) {
		this.id = id;
		this.sequence = sequence;
		this.length = sequence.length;
	}
	
	public Sequence(String id, String sequence) {
		this.id = id;
		this.sequence = sequence.toCharArray();
		this.length = sequence.length();
	}
	
	@Override
	public String getID() {
		return id;
	}
	
	public char[] getSequence() {
		return sequence;
	}
	
	public String getSequenceAsString() {
		return sequence.toString();
	}
	
	@Override
	public int length() {
		return length;
	}

	@Override
	public Character getComp(int i) {
		return sequence[i];
	}
	
}
