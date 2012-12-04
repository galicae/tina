package highscorealignments;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;

public class CathScopList {
	public static ArrayList<CathScopEntry> read(String filepath){
		ArrayList<CathScopEntry> cathscopinfos = new ArrayList<CathScopEntry>();

		try{
			BufferedReader in = new BufferedReader(new FileReader(filepath));
			String line;
			while((line = in.readLine()) != null){
				String[] temp = line.split("\\s");
				String[] cathinfos = temp[temp.length-2].split("\\.");
				String[] scopinfos = temp[temp.length-1].split("\\.");
				cathscopinfos.add(new CathScopEntry(temp[0],cathinfos[0].charAt(0),Integer.parseInt(cathinfos[3]),Integer.parseInt(cathinfos[2]),Integer.parseInt(cathinfos[1]),
						Integer.parseInt(scopinfos[0]),Integer.parseInt(scopinfos[3]),Integer.parseInt(scopinfos[2]),Integer.parseInt(scopinfos[1])));
			}
			in.close();
		} catch(IOException e){
			System.out.println("No pairs input!");
		}
		return cathscopinfos;
	}
}
