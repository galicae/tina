package bioinfo.alignment;

public interface Alignment {
	
	int[][] calcMap();
	int length();
	int getScore();
	String toString();
	Alignable getComponent(int number);
	
}
