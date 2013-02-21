package util;

public class Math {

	public static boolean equals(double a, double b) {
		return (a < b + 0.00001 && a > b - 0.00001);
	}

	public static double dist3D(double x, double y, double z) {
		return java.lang.Math.sqrt(x * x + y * y + z * z);
	}

	public static double zscore(double var, double[] vars) {
		double mean = mean(vars);
		double sd = sd(var, vars, mean);
		return (var - mean) / sd;
	}

	private static double sd(double var, double[] vars, double mean) {
		double sum = 0;
		for (int i = 0; i < vars.length; i++) {
			sum += java.lang.Math.pow(vars[i], 2);
		}
		sum = sum / vars.length;
		return java.lang.Math.sqrt(sum - java.lang.Math.pow(mean, 2));
	}

	private static double mean(double[] var) {
		double sum = 0;
		for (int i = 0; i < var.length; i++) {
			sum += var[i];
		}
		return sum / var.length;
	}
}