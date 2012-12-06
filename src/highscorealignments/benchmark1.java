package highscorealignments;

public class benchmark1 {

	public static void main(String[] args) {
		String[] mode = {"local"};
		String[] matrices = {"dayhoff.mat","BlakeCohenMatrix.mat","THREADERSimilarityMatrix.mat"};
		int[] gapopen = {-9,-12,-15};
		
		for(String mod : mode){
			for(String matrix : matrices){
				for(int go : gapopen){
					BMCathAndScop bm = new BMCathAndScop(args[0], args[1], args[2], args[3]+"/"+matrix, go, -1, mod, mod+"_"+matrix+"_"+go+".bm");
				}
			}
		}
	}
}
