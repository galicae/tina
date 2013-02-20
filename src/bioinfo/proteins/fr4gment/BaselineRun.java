package bioinfo.proteins.fr4gment;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.util.LinkedList;

import bioinfo.alignment.SequenceAlignment;
import bioinfo.alignment.SequenceAlignmentFileReader;
import bioinfo.proteins.PDBEntry;
import bioinfo.proteins.PDBFileReader;
import bioinfo.proteins.fragm3nt.ProteinFragment;
import bioinfo.superpos.Kabsch;
import bioinfo.superpos.PDBReduce;
import bioinfo.superpos.Transformation;


/**
 * class to run a massive test of the baseline method
 * 
 * @author galicae
 * 
 */
public class BaselineRun {

	public static void main(String[] args) throws Exception {
		BufferedWriter wr = new BufferedWriter(new FileWriter("baselineSameFold"));
		SequenceAlignmentFileReader alireader = new SequenceAlignmentFileReader();
		alireader.initSequentialRead("structAlignments/sameFold.scores");
		SequenceAlignment input = null;
		while((input = alireader.nextAlignment()) != null) {
			wr.write(input.toString() + " ");
			LoopBaseline luup = new LoopBaseline(input, "msp_cluster05",
					"cathscop.ids");
			ProteinFragment rr = luup.makePrediction();
			
			PDBFileReader re = new PDBFileReader("/home/galicae/Desktop/STRUCTURES/");
			PDBEntry natStr = re.readFromFolderById("1c3gA01");
			
			ProteinFragment control = new ProteinFragment("control", PDBReduce.reduceSinglePDB(natStr), 4);
			control.setSequence(natStr.getSequenceAsString());
			LinkedList<int[]> cores = new LinkedList<int[]>();
			cores.addAll(luup.getQuCores());
			
			double[][] contr = control.getAllResidues();
			double[][] coord = rr.getAllResidues();
			
			int length = 0;
			for(int[] core: cores) {
				length += core[1] - core[0] + 1;
			}
			
			double[][][] kabschFood = new double[2][length][3];
			
			int c = 0;
			for(int[] core: cores) {
				for(int i = core[0]; i <= Math.min(core[1], contr.length - 1); i++) {
					kabschFood[0][c] = coord[i];
					kabschFood[1][c] = contr[i];
					c++;
				}
			}
			
			Transformation t = Kabsch.calculateTransformation(kabschFood);
			System.out.println(Double.toString(t.getRmsd()));
			wr.write(Double.toString(t.getRmsd()) + "\n");
		}
	}
}
