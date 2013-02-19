package bioinfo.proteins.fr4gment;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.util.LinkedList;

/**
 * this run class is meant to find all distinct pairs in cathscop.inpairs and
 * cathscop.outpairs and write them along with their fold.
 * 
 * @author galicae
 * 
 */
public class MSCatalogue {
	public static void main(String[] args) throws Exception {
		BufferedReader reader = new BufferedReader(new FileReader("cathscop.inpairs"));
		LinkedList<String> alll = new LinkedList<String>();
		LinkedList<String> result = new LinkedList<String>();
		String line = "";
		String[] curr = new String[3];

		while ((line = reader.readLine()) != null) {
			curr = line.split(" ");
			alll.add(curr[0] + " " + curr[6]);
			alll.add(curr[1] + " " + curr[7]);
		}
		reader.close();
		System.out.println("finished inpairs");
//		reader  = new BufferedReader(new FileReader("cathscop.outpairs"));
//		while ((line = reader.readLine()) != null) {
//			curr = line.split(" ");
//			alll.add(curr[0] + " " + curr[6]);
//			alll.add(curr[1] + " " + curr[7]);
//		}
//		
		for (int i = 0; i < alll.size(); i++) {
			System.out.println(i + "/" + alll.size());
			if (!result.contains(alll.get(i)))
				result.add(alll.get(i));
		}
//		reader.close();
		System.out.println("finished outpairs... writing");
		BufferedWriter wr = new BufferedWriter(new FileWriter("cathscop.ids"));
		for(String s: result)
			wr.write(s + "\n");
		
		wr.close();
		System.out.println(result.size());
	}
}
