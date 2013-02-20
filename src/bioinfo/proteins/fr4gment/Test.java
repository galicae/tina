package bioinfo.proteins.fr4gment;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.util.LinkedList;

import bioinfo.Sequence;
import bioinfo.alignment.SequenceAlignment;
import bioinfo.proteins.PDBEntry;
import bioinfo.proteins.PDBFileReader;
import bioinfo.proteins.fragm3nt.ProteinFragment;
import bioinfo.superpos.Kabsch;
import bioinfo.superpos.PDBReduce;
import bioinfo.superpos.Transformation;

/**
 * demonstration of principle of BaselineLoop
 * 
 * @author galicae
 * 
 */
public class Test {
	public static void main(String[] args) throws Exception {
		// System.out.println(RunHelper.execToString("grep -i \"hello world\" ./"));

		// grep 1nltA01
		Sequence seq1 = new Sequence(
				"1nltA01",
				"PQRGKDIKHEISASLEELYKGRTAKLALNKQILVENERKILEVHVEPGMKDGQRIVFKGEADQAPDVIPGDVVFIVSERP");
		Sequence seq2 = new Sequence(
				"1c3gA01",
				"ETVQVNLPVSLEDLFVGKKKSFKIGRKGPHGASEKTQIDIQLKPGWKAGTKITYKNQGDYNPQTGRRKTLQFVIQEKS");

		SequenceAlignment input = new SequenceAlignment(
				seq1,
				seq2,
				"--PQRGKDIKHEIS-----------ASLEELYKGRTAKLALNKQILVENER-----------KILEV-----------------------------HVEPGMKDGQRIVFKGEADQAPDVIPGDVVFIVS----ERP",
				"ET------VQVNLPVSLEDLFVGKK----------------------KSFKIGRKGPHGASEKTQIDIQLKPGWKAGTKITYKNQGDYNPQTGRRK----------------------------TLQFVIQEKS---",
				0.208);

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
		
		double[][][] kabschFood = new double[2][21][3];
		
		int c = 0;
		for(int[] core: cores) {
			for(int i = core[0]; i <= core[1]; i++) {
				kabschFood[0][c] = coord[i];
				kabschFood[1][c] = contr[i];
				c++;
			}
		}
		
		Transformation t = Kabsch.calculateTransformation(kabschFood);
		System.out.println(t.getRmsd());
//		contr = t.transform(contr);
//		control.setCoordinates(contr);
//		
//		try {
//			BufferedWriter wr = new BufferedWriter(
//					new FileWriter("predict.pdb"));
//			wr.write("MODEL        1\n");
//			wr.write(rr.toString());
//			wr.write("ENDMDL\n");
//			wr.write("MODEL        2\n");
//			wr.write(control.toString());
//			wr.write("ENDMDL");
//			wr.close();
//		} catch (Exception e) {
//			e.printStackTrace();
//		}
	}

}
