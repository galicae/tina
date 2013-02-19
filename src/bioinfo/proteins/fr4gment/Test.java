package bioinfo.proteins.fr4gment;

import java.io.BufferedWriter;
import java.io.FileWriter;

import bioinfo.Sequence;
import bioinfo.alignment.SequenceAlignment;
import bioinfo.proteins.fragm3nt.ProteinFragment;
import bioinfo.proteins.fragm3nt.run.RunHelper;

public class Test {
	public static void main(String[] args) throws Exception {
//		System.out.println(RunHelper.execToString("grep -i \"hello world\" ./"));
		
//		 grep 1nltA01 
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

		LoopBaseline luup = new LoopBaseline(input, "msp_cluster05", "cathscop.ids");
		ProteinFragment rr = luup.makePrediction();
		
		try {
			BufferedWriter wr = new BufferedWriter(new FileWriter("predict.pdb"));
			wr.write(rr.toString());
			wr.close();
		}
		catch(Exception e) {
			e.printStackTrace();
		}
	}

	

	
}
