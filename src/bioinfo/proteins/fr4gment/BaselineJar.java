package bioinfo.proteins.fr4gment;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.util.LinkedList;

import joptsimple.OptionParser;
import joptsimple.OptionSet;

import bioinfo.alignment.SequenceAlignment;
import bioinfo.alignment.SequenceAlignmentFileReader;
import bioinfo.proteins.PDBEntry;
import bioinfo.proteins.PDBFileReader;
import bioinfo.proteins.fragm3nt.ProteinFragment;
import bioinfo.superpos.Kabsch;
import bioinfo.superpos.PDBReduce;
import bioinfo.superpos.Transformation;

/**
 * class to run a massive test of the baseline method against MODELLER
 * 
 * @author galicae
 * 
 */
public class BaselineJar {

	// baselineSameFold structAlignments/sameFold.scores msp_cluster05
	// cathscop.ids /home/galicae/Desktop/STRUCTURES/
	public static void main(String[] args) throws Exception {

		String structuralALi_dir = "";
		String multipleSuperpos_dir = "";
		String ids_file = "";
		String protein_files = "";
		String filter = "";
		String output = "";

		OptionParser parser = new OptionParser();

		parser.accepts("ali", "core alignment file").withRequiredArg()
				.ofType(String.class);
		parser.accepts("mss",
				"multiple structural superposition cluster directory")
				.withRequiredArg().ofType(String.class);
		parser.accepts("o", "name of the output file").withRequiredArg()
				.ofType(String.class);
		parser.accepts("pdb", "PDB files directory").withRequiredArg()
				.ofType(String.class);
		parser.accepts("filter", "filter superposition").withRequiredArg()
				.ofType(String.class);
		parser.accepts("ids", "ids file").withRequiredArg()
				.ofType(String.class);
		parser.accepts("h", "show help");

		OptionSet options = parser.parse(args);

		if (options.has("ali")) {
			structuralALi_dir = (String) options.valueOf("ali");
		}
		if (options.has("mss")) {
			multipleSuperpos_dir = (String) options.valueOf("mss");
		}
		if (options.has("o")) {
			output = (String) options.valueOf("o");
		}
		if (options.has("pdb")) {
			protein_files = (String) options.valueOf("pdb");
		}
		if (options.has("filter")) {
			filter = (String) options.valueOf("filter");
		}
		if (options.has("ids")) {
			ids_file = (String) options.valueOf("ids");
		}
		if (options.has("h")) {
			System.out
					.println("========================BASELINE========================\n"
							+ "Baseline loop prediction software based on a multiple \n"
							+ "assembly. Developement abandoned due to poor performance\n"
							+ "and time pressure. Welcome to the help screen.");
			parser.printHelpOn(System.out);
			System.exit(0);
		}

		SequenceAlignmentFileReader alireader = new SequenceAlignmentFileReader();
		alireader.initSequentialRead(structuralALi_dir);
		SequenceAlignment input = null;

		input = alireader.nextAlignment();
		try {
			LoopBaseline luup = new LoopBaseline(input, multipleSuperpos_dir,
					ids_file, filter);
			ProteinFragment rr = luup.makePrediction();

			PDBFileReader re = new PDBFileReader(protein_files);
			PDBEntry natStr = re.readFromFolderById(input.getComponent(1)
					.getId());

			double[][] contr = PDBReduce.reduceSinglePDB(natStr);
			ProteinFragment control = new ProteinFragment("control", contr, 4);
			control.setSequence(natStr.getSequenceAsString());
			LinkedList<int[]> cores = new LinkedList<int[]>();
			cores.addAll(luup.getQuCores());

			double[][] coord = rr.getAllResidues();

			int length = 0;
			for (int[] core : cores) {
				length += core[1] - core[0] + 1;
			}

			double[][][] kabschFood = new double[2][length][3];

			int c = 0;
			for (int[] core : cores) {
				for (int i = core[0]; i <= Math.min(core[1], contr.length - 1); i++) {
					kabschFood[0][c] = coord[i];
					kabschFood[1][c] = contr[i];
					c++;
				}
			}
			
			checkProgress(cores, rr);
			Transformation t = Kabsch.calculateTransformation(kabschFood);
			control.setCoordinates(t.transform(contr));
			contr = control.getAllResidues();
			
			BufferedWriter test = new BufferedWriter(new FileWriter(output));
			test.write("MODEL        1\n");
			test.write(control.toString());
			test.write("ENDMDL\n");
			test.write("MODEL        2\n");
			test.write(rr.toString());
			test.write("ENDMDL\n");
			test.close();
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}
