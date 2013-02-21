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
	public String getId() {
		return id;
	}
	
	public char[] getSequence() {
		return sequence;
	}
	
	public String getSequenceAsString() {
		String result = "";
		for (int i = 0; i < sequence.length; i++) {
			result += sequence[i];
		}
		return result;
	}
	
	@Override
	public int length() {
		return length;
	}

	@Override
	public Character getComp(int i) {
		return sequence[i];
	}
	
	public String toString(){
		String seq = "";
		for(char x: sequence){
			seq += x;
		}
		return seq;
	}
	
	/**
	 * 
	 * @return sequence along with its identifier
	 */
	public String toStringVerbose(){
		String seq = "";
		for(char x: sequence){
			seq += x;
		}
		return id+": "+seq;
	}
	
}
