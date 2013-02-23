package bioinfo.proteins.corecluster;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.HashMap;

public class CoreDefinition {

	public static HashMap<Character, Character> dsspTothreeState = new HashMap<Character, Character>();

	public static void init() {
		dsspTothreeState = new HashMap<Character, Character>();
		dsspTothreeState.put('H', 'H');
		dsspTothreeState.put('G', 'H');
		dsspTothreeState.put('I', 'H');
		dsspTothreeState.put('B', 'E');
		dsspTothreeState.put('E', 'E');
		dsspTothreeState.put('T', 'C');
		dsspTothreeState.put('C', 'C');
		dsspTothreeState.put('S', 'C');
		dsspTothreeState.put('-', '-');
	}

	public static HashMap<String, char[]> parseDsspToThreeState(String f) {
		init();
		HashMap<String, char[]> dssp = new HashMap<String, char[]>();
		BufferedReader bf;
		try {
			bf = new BufferedReader(new FileReader(f));
			String s;

			String name = "";
			while ((s = bf.readLine()) != null) {

				if (s.startsWith(">")) {
					name = s.substring(1);
				} else {
					char[] dsspchar = s.toCharArray();
					char[] threestate = new char[s.length() + 1];
					for (int i = 0; i < dsspchar.length; i++) {
						threestate[i] = dsspTothreeState.get(dsspchar[i]);
					}
					threestate[threestate.length - 1] = 'C';
					dssp.put(name, threestate);
				}

			}
			bf.close();
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		return dssp;
	}
}
