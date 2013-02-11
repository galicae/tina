package test;

import java.io.BufferedReader;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.Arrays;

public class ReprofTest {
	
	public static void main(String[] args) throws Exception {
		reprof("lib/reprof", "~/Desktop/example.fasta", "/home/galicae/Desktop");
		
	}
	
	public static void reprof(String reprofDir, String input, String output) throws Exception {
		if(reprofDir.endsWith("/"))
			reprofDir = reprofDir.substring(0, reprofDir.length()-1);
		
		String model = " --model " + reprofDir + "/share/";
		input = " --input " + input;
//		output = " --output " + output;
//		input = "";
		output = "";
		model = "";
		
		String call = reprofDir + "/bin/reprof" + input + output + model;
		System.out.println(call);
		Process process = new ProcessBuilder(call).start();
		InputStream is = process.getInputStream();
		InputStreamReader isr = new InputStreamReader(is);
		BufferedReader br = new BufferedReader(isr);
		String line;

		while ((line = br.readLine()) != null) {
		  System.out.println(line);
		}

	}
}
