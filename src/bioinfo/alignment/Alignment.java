package bioinfo.alignment;

public interface Alignment {
	
	int[][] getAlignedResidues();
	int[][] calcMap();
	int length();
	double getScore();
	String toString();
	String toStringVerbose();
	Alignable getComponent(int number);
	Alignment duplicate();
	
}
