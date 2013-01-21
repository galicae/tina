package test;

import java.io.BufferedReader;
import java.io.FileReader;

import bioinfo.pdb.PDBFile;
import bioinfo.proteins.PDBEntry;

public class toStringTest {

	public static void main(String[] args) {
		try {
			BufferedReader br = new BufferedReader(new FileReader("prot_list"));
			String line;
			while((line = br.readLine()) != null) {
				PDBFile.downloadPDB(line, "./proteins/");
				System.out.println("got " + line);
			}
			br.close();
		}
		catch(Exception e) {
			e.printStackTrace();
		}
		
//		PDBFile.downloadPDB("1BMV", "./proteins/");
	}
}
