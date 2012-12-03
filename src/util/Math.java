package util;

public class Math {

	public static boolean equals(double a, double b){
		return (a < b +0.00001 && a > b - 0.00001);
	}
	
	public static double dist3D(double x, double y, double z){
		return java.lang.Math.sqrt(x*x+y*y+z*z);
	}
}
