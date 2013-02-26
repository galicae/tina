package bioinfo.proteins.fragm3nt.run;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.LinkedList;

import joptsimple.OptionParser;
import joptsimple.OptionSet;
import bioinfo.proteins.PDBEntry;
import bioinfo.proteins.fragm3nt.ClusterAnalysis;
import bioinfo.proteins.fragm3nt.FragmentCluster;
import bioinfo.proteins.fragm3nt.ProteinFragment;
import bioinfo.proteins.fragm3nt.assembler.Assembler;

/**
 * run class to be turned into a jar
 * 
 * @author galicae
 * 
 */
public class Fragm3ntRun {

	public static void main(String[] args) throws IOException {
		LinkedList<FragmentCluster> clusters = new LinkedList<FragmentCluster>();
		String output = "./test";
		int fragLength = 8;
		int extent = 5;
		String query = "input";

		OptionParser parser = new OptionParser();

		parser.accepts("s", "sequence to be predicted").withRequiredArg()
				.ofType(String.class);
		parser.accepts("c", "cluster folder").withRequiredArg()
				.ofType(String.class);
		parser.accepts("o", "name of the output file").withRequiredArg()
				.ofType(String.class);
		parser.accepts("e",
				"extent to which fragments must overlap (a minimum of 4 is recommended)")
				.withOptionalArg().ofType(Integer.class);
		parser.accepts("h", "show help");

		OptionSet options = parser.parse(args);

		if (options.has("h")) {
			System.out
					.println("========================FRAGM3NT========================\n"
							+ "Protein backbone prediction software based on fragment \n"
							+ "assembly. Developement abandoned due to poor performance\n"
							+ "and time pressure. Welcome to the help screen.");
			parser.printHelpOn(System.out);
			System.exit(0);
		}

		if (options.has("s")) {
			query = (String) options.valueOf("s");
		}
		if (options.has("c")) {
			ClusterAnalysis c = new ClusterAnalysis(
					(String) options.valueOf("c"));
			clusters = c.getClusters();
		}
		if(options.has("o")) {
			output = (String) options.valueOf("o");
		}
		if(options.has("e")) {
			extent = (Integer) options.valueOf("e");
			if(extent < 1) {
				System.out.println("An invalid extent was detected. Defaulting to 4.");
				extent = 4;
			}
		}

		Assembler ass = new Assembler(fragLength);
		ProteinFragment result = ass.predictStructure(query, clusters, extent);

		PDBEntry resultEntry = result.toPDB();

		try {
			BufferedWriter wr = new BufferedWriter(new FileWriter(output));
			wr.write(resultEntry.getAtomSectionAsString());
			wr.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
}
