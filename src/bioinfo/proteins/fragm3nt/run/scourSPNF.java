package bioinfo.proteins.fragm3nt.run;

import java.io.BufferedReader;
import java.io.FileReader;
import java.util.LinkedList;

/**
 * the purpose of this class is to find the number of distinct protein ids in
 * file spnfTM, which contains the pairs from cathscop.inpairs that belong to
 * the same family but not the same superfamily and have a TMalign score upwards
 * of 0.5
 * 
 * @author galicae
 * 
 */
public class scourSPNF {

	public static void main(String[] args) throws Exception {
		BufferedReader reader = new BufferedReader(new FileReader("spnfTM"));
		LinkedList<String> alll = new LinkedList<String>();
		LinkedList<String> result = new LinkedList<String>();
		String line = "";
		String[] curr = new String[3];

		while ((line = reader.readLine()) != null) {
			curr = line.split(" ");
			alll.add(curr[0]);
			alll.add(curr[1]);
		}

		reader.close();

		for (int i = 0; i < alll.size(); i++) {
			if (!result.contains(alll.get(i)))
				result.add(alll.get(i));
		}

		System.out.println(result.size());
	}
}
