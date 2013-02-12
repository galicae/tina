package bioinfo.proteins;

import java.io.BufferedReader;
import java.io.FileReader;
import java.util.ArrayList;

import org.apache.commons.exec.CommandLine;
import org.apache.commons.exec.DefaultExecutor;

public class ReProf {

	/**
	 * this function calls the reprof binary from its source destination (it is
	 * imperative that the structure of the folder is not compromised)
	 * 
	 * @param reprofDir
	 *            the folder reprof, containing the bin folder with the binary
	 * @param input
	 *            the input fasta file
	 * @param out
	 *            the output directory
	 * @throws Exception
	 */
	public static void predictSecStruct(String reprofDir, String input,
			String out) throws Exception {
		if (reprofDir.endsWith("/"))
			reprofDir = reprofDir.substring(0, reprofDir.length() - 1);

		String model = reprofDir + "/reprof/share";
		// input = "--input" + input;

		String call = reprofDir + "/bin/reprof";
		call = call + " --input " + input + " --model " + model + " --single "
				+ " --out " + out;
		// , input, output, model, " --single"
		CommandLine cmdLine = CommandLine.parse(call);
		DefaultExecutor executor = new DefaultExecutor();
		executor.execute(cmdLine);
	}

	/**
	 * 
	 * this function calls the reprof binary from its source destination (it is
	 * imperative that the structure of the folder is not compromised)
	 * 
	 * @param reprofDir
	 *            the folder reprof, containing the bin folder with the binary
	 * @param input
	 *            the input fasta file
	 * @throws Exception
	 * 
	 * @return an array of SecStructThree
	 */
	public static SecStructThree[] predictSecStructVerb(String reprofDir,
			String input, String out) {
		try {
			if (reprofDir.endsWith("/"))
				reprofDir = reprofDir.substring(0, reprofDir.length() - 1);

			String model = reprofDir + "/reprof/share";
			// input = "--input" + input;

			String call = reprofDir + "/bin/reprof";
			call = call + " --input " + input + " --model " + model
					+ " --single " + " --out " + out;
			// , input, output, model, " --single"
			CommandLine cmdLine = CommandLine.parse(call);
			DefaultExecutor executor = new DefaultExecutor();
			executor.execute(cmdLine);

			BufferedReader br = new BufferedReader(new FileReader(out + "_ORI"));
			String line;
			ArrayList<SecStructThree> resList = new ArrayList<SecStructThree>();
			while ((line = br.readLine()) != null) {
				if(Character.isDigit(line.charAt(0))){
					String[] lineArr = line.split("\t");
					resList.add();
				}
			}
		} catch (Exception e) {

		}
		return null;

	}
}
