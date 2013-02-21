package bioinfo.proteins;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.util.ArrayList;

import org.apache.commons.exec.CommandLine;
import org.apache.commons.exec.DefaultExecutor;

import bioinfo.alignment.kerbsch.temp.SecStructScores;

public class ReProf {

	public static SecStructThree[] predictSecStruct(PDBEntry pdb, String path,
			String reprofDir) {
		String input = createFasta(pdb, path);
		return predictSecStructVerb(reprofDir, input, path);
	}

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
	private static void predictSecStruct(String reprofDir, String input,
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
	private static SecStructThree[] predictSecStructVerb(String reprofDir,
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
				if (Character.isDigit(line.charAt(0))) {
					String[] lineArr = line.split("\t");
					resList.add(SecStructThree.defSecStructThree(lineArr[2]
							.charAt(0)));
				}
			}
			br.close();
			SecStructThree[] result = new SecStructThree[resList.size()];
			for (int i = 0; i < resList.size(); i++) {
				result[i] = resList.get(i);
			}
			File output = new File(out + "_ORI");
			output.delete();
			File fasta = new File(input);
			fasta.delete();
			return result;

		} catch (Exception e) {

		}
		return null;
	}

	/**
	 * this function creates a FASTA file from a PDBEntry and returns its
	 * location
	 * 
	 * @param pdb
	 *            the PDBEntry to parse
	 * @param path
	 *            the path where the file is to be placed
	 * @return the path to the FASTA file
	 */
	private static String createFasta(PDBEntry pdb, String path) {
		StringBuilder sb = new StringBuilder();
		for (int i = 0; i < pdb.length(); i++) {
			sb.append(pdb.getAminoAcid(i).getName());
			if (i % 80 == 0 && i > 0)
				sb.append("\n");
		}
		sb.insert(0, ">" + pdb.getId() + "\n");

		try {
			BufferedWriter wr = new BufferedWriter(new FileWriter(path
					+ pdb.getId() + ".fasta"));
			wr.write(sb.toString());
			wr.close();
		} catch (Exception e) {
			e.printStackTrace();
		}
		return path + pdb.getId() + ".fasta";
	}
}
