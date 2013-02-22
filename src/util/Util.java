package util;

import static java.lang.Math.sqrt;

import java.util.Locale;

import bioinfo.proteins.AminoAcid;
import bioinfo.proteins.Atom;
import bioinfo.proteins.AtomType;


/**
 * The class Util contains some useful methods
 * @author gobi_12_4
 * @lastchange 2013-02-20
 */
public class Util {
	
	public static void printDoubleArray(double[][] array, int prec){
		for(int i = 0; i != array.length; i++){
			for(int j = 0; j != array[0].length; j++){
				System.out.print(String.format(Locale.US,"%."+prec+"f", array[i][j])+"\t");
			}
			System.out.println();
		}
	}
	
	public static void printDoubleArray(double[] array, int prec){
		for(int i = 0; i != array.length; i++){
			System.out.print(String.format(Locale.US,"%."+prec+"f", array[i])+"\t");
		}
		System.out.println();
	}
	
	public static void printIntegerArray(int[][] array){
		for(int i = 0; i != array.length; i++){
			for(int j = 0; j != array[0].length; j++){
				System.out.print(array[i][j]+"\t");
			}
			System.out.println();
		}
	}
	
	public static void printIntegerArray(int[] array){
		for(int i = 0; i != array.length; i++){
			System.out.print(array[i]+"\t");
		}
		System.out.println();
	}
	
	/**
	 * makes sure that the given folder String ends with a "/"
	 * @param folder
	 * @return
	 */
	public static String checkedFolder(String folder){
		if(folder.endsWith("/")){
			return folder;
		} else{
			return folder+"/";
		}
	}
	
	/**
	 * flips a char[] on itself
	 * runtime: n. memory: 2n
	 * @param in the character array in question
	 * @return the reversed array
	 */
	public static char[] flip(char[] in) {
		char[] out = new char[in.length];
		for (int i = in.length - 1; i >= 0; i--) {
			out[out.length - 1 - i] = in[i];
		}
		return out;
	}
	
	/**
	 * calculates the euklidian distance between two AminoAcids
	 * 
	 * @param a
	 *            an AmoniAcid
	 * @param b
	 *            another AminoAcid
	 * @return the euklidian distance between the two AminoAcid's C alpha atoms
	 */
	public static double calcDistance(AminoAcid a, AminoAcid b) {
		Atom caa = a.getAtomByType(AtomType.CA);
		Atom cab = b.getAtomByType(AtomType.CA);
		if (caa != null && cab != null) {
			return calcDistance(caa, cab);
		}
		return Double.NaN;
	}

	/**
	 * calculates the distance between two atoms
	 * 
	 * @param a
	 *            an Atom
	 * @param b
	 *            another Atom
	 * @return the euklidian distance between two Atoms
	 */
	public static double calcDistance(Atom a, Atom b) {
		double[] apos = a.getPosition();
		double[] bpos = b.getPosition();
		double[] dis = { apos[0] - bpos[0], apos[1] - bpos[1],
				apos[2] - bpos[2] };
		return sqrt((dis[0] * dis[0]) + (dis[1] * dis[1])
				+ (dis[2] * dis[2]));
	}
	
}
