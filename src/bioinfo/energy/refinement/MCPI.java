package bioinfo.energy.refinement;

import java.util.Random;

public class MCPI {

	public static void main(String[] args) {
		int size = Integer.parseInt(args[0]);
		double inside = 0.0d;
		Random rand = new Random();
		double x;
		double y;
		double d;
		for(int i = 0; i != size; i++){
			x = rand.nextDouble();
			y = rand.nextDouble();
			d = Math.sqrt((x*x)+(y*y));
			if(d <=  1.0d){
				inside++;
			}
		}
		System.out.println(inside);
		System.out.println(4*(inside/size));
	}

}
