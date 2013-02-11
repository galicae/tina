package test;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.util.LinkedList;

public class WaitBar {
	public static void main(String[] args) {
		LinkedList<String> idList = new LinkedList<String>();
		int count = 0;
		try {
			BufferedReader reader = new BufferedReader(new FileReader(
					"list.txt"));
			String line = reader.readLine();
			String[] lArr = new String[2];
			while (line != null) {
				lArr = line.split(" ");
				if (lArr[1].equals("b.1.1.4")){
					count++;
					idList.add(lArr[0]);
				}
				line = reader.readLine();
			}
			reader.close();
			System.out.println("finished reading  " + count);
		} catch (Exception e) {
			e.printStackTrace();
		}

		try {
			BufferedWriter wr = new BufferedWriter(new FileWriter("idList.txt"));
			for (int i = 0; i < idList.size(); i++) {
				wr.write(idList.get(i) + "\n");
			}
			wr.close();
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}