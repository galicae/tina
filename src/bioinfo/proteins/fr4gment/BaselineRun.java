package bioinfo.proteins.fr4gment;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.util.LinkedList;

import bioinfo.alignment.SequenceAlignment;
import bioinfo.alignment.SequenceAlignmentFileReader;
import bioinfo.proteins.PDBEntry;
import bioinfo.proteins.PDBFileReader;
import bioinfo.proteins.fragm3nt.ProteinFragment;
import bioinfo.proteins.fragm3nt.run.RunHelper;
import bioinfo.superpos.Kabsch;
import bioinfo.superpos.PDBReduce;
import bioinfo.superpos.Transformation;

/**
 * class to run a massive test of the baseline method against MODELLER
 * 
 * @author galicae
 * 
 */
public class BaselineRun {
	// result_filename structuralALi_dir multipleSuperpos_dir ids_file
	// protein_files
	// baselineSameFold structAlignments/sameFold.scores msp_cluster05
	// cathscop.ids /home/galicae/Desktop/STRUCTURES/
	public static void main(String[] args) throws Exception {
		BufferedWriter wr = new BufferedWriter(new FileWriter(args[0]));
		System.out.println("seque_1 seque_2 bsRMSD mdRMSD lL #c");
		wr.write("seque_1 seque_2 bsRMSD mdRMSD lL #c\n");
		SequenceAlignmentFileReader alireader = new SequenceAlignmentFileReader();
		alireader.initSequentialRead(args[1]);
		SequenceAlignment input = null;

		while ((input = alireader.nextAlignment()) != null) {
			try {
				// if (counter > 0) {
				// counter--;
				// continue;
				// }

				ProteinFragment modellerPred = runModeller(input, args[5]);

				LoopBaseline luup = new LoopBaseline(input, args[2], args[3], args[6]);
				ProteinFragment rr = luup.makePrediction();

				PDBFileReader re = new PDBFileReader(args[4]);
				PDBEntry natStr = re.readFromFolderById(input.getComponent(1)
						.getId());

				double[][] contr = PDBReduce.reduceSinglePDB(natStr);
				ProteinFragment control = new ProteinFragment("control", contr,
						4);
				control.setSequence(natStr.getSequenceAsString());
				LinkedList<int[]> cores = new LinkedList<int[]>();
				cores.addAll(luup.getQuCores());

				double[][] coord = rr.getAllResidues();
				double[][] mod = modellerPred.getAllResidues();

				int length = 0;
				for (int[] core : cores) {
					length += core[1] - core[0] + 1;
				}

				double[][][] kabschFood = new double[2][length][3];

				int c = 0;
				for (int[] core : cores) {
					for (int i = core[0]; i <= Math.min(core[1],
							contr.length - 1); i++) {
						kabschFood[0][c] = coord[i];
						kabschFood[1][c] = contr[i];
						c++;
					}
				}

				Transformation t = Kabsch.calculateTransformation(kabschFood);
				control.setCoordinates(t.transform(contr));
				contr = control.getAllResidues();

				// BufferedWriter test = new BufferedWriter(new
				// FileWriter("test"));
				// test.write("MODEL 1\n");
				// test.write(control.toString());
				// test.write("ENDMDL\n");
				// test.write("MODEL 2\n");
				// test.write(rr.toString());
				// test.write("ENDMDL\n");
				// test.close();

				for (int i = 1; i < cores.size(); i++) {
					int[] prev = cores.get(i - 1);
					int[] next = cores.get(i);

					int loopLength = next[0] - prev[1] + 1;
					double baseRmsd = calcLoopGDT_TS(coord, contr, prev, next);

					if (baseRmsd >= 0) {

						kabschFood = new double[2][c][3];
						c = 0;
						for (int[] core : cores) {
							for (int j = core[0]; j <= Math.min(core[1],
									contr.length - 1); j++) {
								kabschFood[0][c] = contr[j];
								kabschFood[1][c] = mod[j];
								c++;
							}
						}
						Transformation tr = Kabsch
								.calculateTransformation(kabschFood);
						mod = tr.transform(mod);
						modellerPred.setCoordinates(mod);
						double modRmsd = calcLoopGDT_TS(mod, contr, prev, next);
						String seq1 = input.getComponent(0).getId();
						String seq2 = input.getComponent(1).getId();
						wr = new BufferedWriter(new FileWriter(args[0], true));
						System.out.printf("%s %s %.4f %.4f %d %d\n", seq1,
								seq2, baseRmsd, modRmsd, loopLength,
								cores.size());
						wr.write(String.format("%s %s %.4f %.4f %d %d\n", seq1,
								seq2, baseRmsd, modRmsd, loopLength,
								cores.size()));
						wr.close();
					}
				}
				BufferedWriter test = new BufferedWriter(new FileWriter(
						"/home/p/papadopoulos/" + args[5] + "/"
								+ input.getComponent(0).getId() + "_"
								+ input.getComponent(1).getId() + ".pdb"));
				test.write("MODEL 1\n");
				test.write(control.toString());
				test.write("ENDMDL\n");
				test.write("MODEL 2\n");
				test.write(rr.toString());
				test.write("ENDMDL\n");
				test.write("MODEL 3\n");
				test.write(modellerPred.toString());
				test.write("ENDMDL\n");
				test.close();
			} catch (Exception e) {
				e.printStackTrace();
				continue;
			}
		}
		System.out.println("finished with alignments");

	}

	/**
	 * checks if a certain loop region has an rmsd to offer
	 * 
	 * @param input
	 * @param prev
	 * @param next
	 * @return
	 */
	private static boolean checkPrediction(double[][] input, int[] prev,
			int[] next) {
		for (int j = prev[1] + 1; j < next[0]; j++) {
			if (input[j][0] == input[j][1] && input[j][0] == 0)
				continue;
			return true;
		}
		return false;
	}

	/**
	 * calculates the mean rmsd of the aligned regions
	 * 
	 * @param coord
	 * @param contr
	 * @param cores
	 * @return
	 */
	private static double calcLoopGDT_TS(double[][] coord, double[][] contr,
			int[] prev, int[] next) {
		if (!checkPrediction(coord, prev, next))
			return -1;
		double[][] a = new double[next[0] - prev[1] + 1][3];
		double[][] b = new double[next[0] - prev[1] + 1][3];
		for (int i = prev[1]; i <= Math.min(next[0], contr.length - 1); i++) {
			a[i - prev[1]] = coord[i];
			b[i - prev[1]] = contr[i];
		}
		double rmsd = calcGDT_TS(a, b);

		return rmsd;
	}

	private static double calcRmsd(double[][] p, double[][] q) {
		double rmsd = 0;
		double tx = 0, ty = 0, tz = 0;
		for (int i = 0; i < p.length; i++) {
			tx = p[i][0] - q[i][0];
			ty = p[i][1] - q[i][1];
			tz = p[i][2] - q[i][2];

			tx *= tx;
			ty *= ty;
			tz *= tz;

			rmsd = (tx + ty + tz) * (1.0 * p.length);
		}
		rmsd = Math.sqrt(rmsd);
		return rmsd;
	}

	
	private static double calcGDT_TS(double[][] p, double[][] q) {
		double result = 0;
		int counter1 = 0;
		int counter2 = 0;
		int counter4 = 0;
		int counter8 = 0;
		double distance = 0;
		
		for (int i = 0; i < p.length; i++) {
			distance = euclideanDistance(p[i], q[i]);
			if(distance < 1)
				counter1++;
			if(distance < 2)
				counter2++;
			if(distance < 4)
				counter4++;
			if(distance < 8)
				counter8++;
		}
		result = (counter1 + counter2 + counter4 + counter8) / (4.0 * p.length);
		return result;
	}
	
	
	private static double euclideanDistance(double[] x, double[] y) {
		double result = 0;
		for (int i = 0; i < x.length; i++) {
			result += Math.pow(x[i] - y[i], 2);
		}
		return Math.sqrt(result);
	}
	
	/**
	 * this method encapsulates the MODELLER call and reads the prediction PDB
	 * file
	 * 
	 * @param input
	 * @return
	 */
	private static ProteinFragment runModeller(SequenceAlignment input,
			String dir) throws Exception {
		// first write the alignment to a file
		BufferedWriter aliw = new BufferedWriter(new FileWriter(
				"/home/p/papadopoulos/" + dir + "/test.ali"));
		// String newline = System.getProperty("line.separator");

		String[] ali = input.toStringVerbose().split("\n");
		ali[1] = ali[1].replace(" ", "");
		ali[2] = ali[2].replace(" ", "");
		// System.out.println(ali[1]);
		aliw.write(ali[2] + "\n");
		aliw.write(ali[1] + "\n");
		aliw.close();
		String command = "/home/proj/biosoft/PROTEINS/scripts/model_full.sh -i /home/p/papadopoulos/"
				+ dir
				+ "/test.ali -S /home/proj/biosoft/PROTEINS/CATHSCOP/STRUCTURES/ -o /home/p/papadopoulos/"
				+ dir + "/" + input.getComponent(1).getId() + ".pdb";
		RunHelper.execToString(command);

		PDBFileReader read = new PDBFileReader();
		PDBEntry model = read.readPDBFromFile("/home/p/papadopoulos/" + dir
				+ "/" + input.getComponent(1).getId() + ".pdb");
		double[][] modelCoordinates = PDBReduce.reduceSinglePDB(model);
		String sequence = model.getSequenceAsString();
		ProteinFragment result = new ProteinFragment(model.getId(), sequence,
				modelCoordinates, sequence.length());
		return result;
	}
}
