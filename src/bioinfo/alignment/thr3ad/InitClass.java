package bioinfo.alignment.thr3ad;

import java.lang.reflect.Field;

public class InitClass {
	public double[] calcHydropathyMatrix() {
		Field[] vars = HydropathyScores.class.getFields();
		double[] hydropathyScore = new double[26];
		for (Field f : vars) {
			char cur = f.getName().charAt(0);
			try {
				hydropathyScore[cur - 65] = f.getDouble(f);
			} catch (IllegalArgumentException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			} catch (IllegalAccessException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}
		return hydropathyScore;
	}
}
