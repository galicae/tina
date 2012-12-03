package util;

import java.util.Locale;

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
	
	public static String checkedFolder(String folder){
		if(folder.endsWith("/")){
			return folder;
		} else{
			return folder+"/";
		}
	}

	
}
