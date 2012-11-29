package bioinfo.alignment;

public interface Alignable {
	
	/**
	 * 
	 * @param i integer defining position of component in Alignable Sequence
	 * @return component at position i
	 */
	public Object getComp(int i);

	/**
	 * 
	 * @return length of components to align
	 */
	public int length();
	
	/**
	 * 
	 * @return identifier of alignable sequence
	 */
	public String getID();
	
}
