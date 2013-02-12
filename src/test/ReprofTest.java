package test;

import java.io.BufferedReader;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.Arrays;

public class ReprofTest {

	public static void main(String[] args) throws Exception {
		reprof("lib/reprof", "/home/galicae/Desktop/example.fasta",
				"/home/galicae/Desktop");

	}

	public static void reprof(String reprofDir, String input, String output)
			throws Exception {
		if (reprofDir.endsWith("/"))
			reprofDir = reprofDir.substring(0, reprofDir.length() - 1);

		String model = reprofDir + "/reprof/share";
//		input = "--input" + input;
		output = "--output" + output;

		String call = reprofDir + "/bin/reprof";
		// , input, output, model, " --single"
		System.out.println("");
		Process process = new ProcessBuilder(call, "--input", input, "--output", output, "--model", model, "--single").start();
		InputStream is = process.getInputStream();
		InputStreamReader isr = new InputStreamReader(is);
		BufferedReader br = new BufferedReader(isr);
		String line;

		while ((line = br.readLine()) != null) {
			System.out.println(line);
		}

	}
}
