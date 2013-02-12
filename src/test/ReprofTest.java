package test;

import org.apache.commons.exec.CommandLine;
import org.apache.commons.exec.DefaultExecutor;

import bioinfo.proteins.ReProf;
import bioinfo.proteins.SecStructThree;

public class ReprofTest {

	public static void main(String[] args) throws Exception {
		SecStructThree[] d = ReProf.predictSecStructVerb("lib/reprof",
				"/home/galicae/Desktop/example.fasta",
				"/home/galicae/Desktop/example");
		for(int i = 0; i < d.length; i++) {
			System.out.print(d[i].getCharRepres());
		}
	}

	public static void reprof(String reprofDir, String input, String out)
			throws Exception {
		if (reprofDir.endsWith("/"))
			reprofDir = reprofDir.substring(0, reprofDir.length() - 1);

		String model = reprofDir + "/reprof/share";
		// input = "--input" + input;
		// output = "--output" + output;

		String call = reprofDir + "/bin/reprof";
		call = call + " --input " + input + " --model " + model + " --single "
				+ " --out" + " /home/galicae/Desktop/example";
		// , input, output, model, " --single"
		CommandLine cmdLine = CommandLine.parse(call);
		DefaultExecutor executor = new DefaultExecutor();
		int exitValue = executor.execute(cmdLine);
	}
}
