package bioinfo.proteins.fragm3nt.run;

import java.io.BufferedWriter;
import java.io.FileWriter;

import bioinfo.Sequence;
import bioinfo.alignment.SequenceAlignment;
import bioinfo.alignment.gotoh.FreeshiftSequenceGotoh;
import bioinfo.alignment.matrices.QuasarMatrix;
import bioinfo.proteins.PDBEntry;
import bioinfo.proteins.PDBFileReader;
import bioinfo.proteins.fragm3nt.ProteinFragment;
import bioinfo.proteins.fragm3nt.assembler.AlignmentAssembler;
import bioinfo.superpos.Kabsch;
import bioinfo.superpos.PDBReduce;
import bioinfo.superpos.Transformation;

/**
 * simple test class recapping the kabsch pipeline
 * 
 * @author galicae
 * 
 */
public class KabschTest {

	public static void main(String[] args) throws Exception {
		String[] seqIds = { "1chmB01", "1kp0A01" };
		String[] stringSeq = new String[seqIds.length];

		AlignmentAssembler ass = new AlignmentAssembler(8);
		PDBFileReader read = new PDBFileReader(
				"/home/galicae/Desktop/STRUCTURES/");
		PDBEntry pdb = read.readFromFolderById("1chmB01");
		PDBEntry pdb1 = read.readFromFolderById("1kp0A01");

		FreeshiftSequenceGotoh got = new FreeshiftSequenceGotoh(-15, -3,
				QuasarMatrix.DAYHOFF_MATRIX);
		Sequence seq1 = new Sequence(pdb.getID(), pdb.getSequence());
		Sequence seq2 = new Sequence(pdb1.getID(), pdb1.getSequence());
		SequenceAlignment ali = got.align(seq1, seq2);

		BufferedWriter wr2 = new BufferedWriter(new FileWriter("fastaTest.pdb"));

		// double[][][] kabschFood = PDBReduce.reduce(ali, pdb, pdb1);
		// Transformation t = Kabsch.calculateTransformation(kabschFood);
		// System.out.println(t.getRmsd());

		ProteinFragment pdbFrag = new ProteinFragment("real", "D",
				new double[0][0], 8);
		pdbFrag.append(PDBReduce.reduceSinglePDB(pdb), pdb.getSequence());

		ProteinFragment pdb1Frag = new ProteinFragment("real", "D",
				new double[0][0], 8);
		// // pdb1 = t.transform(pdb1);
		pdb1Frag.append(PDBReduce.reduceSinglePDB(pdb1), pdb1.getSequence());

		wr2.write("MODEL        1\n");
		wr2.write(pdbFrag.toString());
		wr2.write("ENDMDL\n");
		wr2.write("MODEL        2\n");
		wr2.write(pdb1Frag.toString());
		wr2.write("ENDMDL\n");
		wr2.close();
	}

}
