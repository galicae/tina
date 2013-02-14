package test;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.util.LinkedList;

public class toStringTest {

	public static void main(String[] args) throws Exception{
		LinkedList<String> result = new LinkedList<String>();
		BufferedReader r1 = new BufferedReader(new FileReader("famTestResults"));
		BufferedReader r2 = new BufferedReader(new FileReader("actVSmaxVSrmsd"));
		String line = r1.readLine();
		while((line = r1.readLine()) != null) {
			String current = line.split("\t")[0];
			String line2 = "";
			while(!(line2 = r2.readLine()).split("\t")[0].equals(current)) {
				
			}
			result.add(line2 + "\t" + line.split("\t")[1] + "\n");
		}
		
		r1.close();
		r2.close();
		
		BufferedWriter w = new BufferedWriter(new FileWriter("result"));
		for(int i = 0; i<result.size(); i++) {
			w.write(result.get(i));
		}
		w.close();
	}
}
