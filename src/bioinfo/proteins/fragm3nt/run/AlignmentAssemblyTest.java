package bioinfo.proteins.fragm3nt.run;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.util.LinkedList;

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

public class AlignmentAssemblyTest {

	public static void main(String[] args) {
		// String[] seqIds = { "1cqxB02", "1gvhA02"};
		String[] seqIds = { "1f0yA02", "3hdhC02" }; // 0.77
		// String[] seqIds = { "1m3uE00", "1o68B00"}; //0.53
		String[] stringSeq = new String[seqIds.length];

		// read sequences (strings)
		try {
			BufferedReader r = new BufferedReader(new FileReader(
					"domains.seqlib"));
			String line = "";
			int count = 0;
			while (count < seqIds.length && (line = r.readLine()) != null) {
				for (int i = 0; i < seqIds.length; i++) {
					if (line.startsWith(seqIds[i])) {
						stringSeq[i] = line.split(":")[1];
						count++;
						continue;
					}
				}
			}
			r.close();
		} catch (Exception e) {
			e.printStackTrace();
		}

		LinkedList<Sequence> seqs = new LinkedList<Sequence>();
		for (int i = 0; i < stringSeq.length; i++) {
			seqs.add(new Sequence(seqIds[i], stringSeq[i]));
		}

		int extent = 4;
		AlignmentAssembler ass = new AlignmentAssembler(8);
		ProteinFragment result = ass.predStrucFromAl(seqs, extent,
				"/home/galicae/Desktop/STRUCTURES/");
		PDBFileReader read = new PDBFileReader(
				"/home/galicae/Desktop/STRUCTURES/");
		PDBEntry pdb = read.readFromFolderById("1f0yA02");
		String query = "";
		for (int i = 0; i < pdb.length(); i++) {
			query += pdb.getAminoAcid(i).getName();
		}
		double[][][] kabschFood = new double[2][][];

		try {
			BufferedWriter wr2 = new BufferedWriter(new FileWriter(
					"fastaTest.pdb"));
			kabschFood = new double[2][result.getAllResidues().length][3];
			kabschFood[0] = result.getAllResidues();
			kabschFood[1] = PDBReduce.reduceSinglePDB(pdb);
			Transformation t = Kabsch.calculateTransformation(kabschFood);
			kabschFood[1] = t.transform(kabschFood[1]);
			System.out.println(t.getRmsd());

			ProteinFragment prot = new ProteinFragment("real", "D",
					new double[0][0], 8);
			prot.setSequence(query);
			prot.append(PDBReduce.reduceSinglePDB(pdb), "");
			pdb = t.transform(pdb);
			prot.setCoordinates(kabschFood[1]);

			wr2.write("MODEL        1\n");
			wr2.write(result.toString(ass.getFragments(), extent));
			wr2.write("ENDMDL\n");
			wr2.write("MODEL        2\n");
			wr2.write(prot.toString());
			wr2.write("ENDMDL\n");
			
//			kabschFood[0] = prot.getAllResidues();
			Sequence seq1 = new Sequence(pdb.getID(), pdb.getSequence());
			for (int i = 1; i < seqs.size(); i++) {
				PDBEntry pdb1 = read.readFromFolderById(seqs.get(i).getID());
				FreeshiftSequenceGotoh got = new FreeshiftSequenceGotoh(-15, -3, QuasarMatrix.DAYHOFF_MATRIX);
				Sequence seq2 = new Sequence(pdb1.getID(), pdb1.getSequence());
				SequenceAlignment ali = got.align(seq1, seq2);
				System.out.println(ali.toStringVerbose());
				
				prot = new ProteinFragment(pdb1.getID(), "D", new double[0][0], 8);
				prot.setSequence(pdb1.getSequence());
				prot.append(PDBReduce.reduceSinglePDB(pdb1), "");
				
				kabschFood[1] = prot.getAllResidues();
				kabschFood = PDBReduce.reduce(ali, pdb, pdb1);
				t = Kabsch.calculateTransformation(kabschFood);
				kabschFood[1] = t.transform(kabschFood[1]);
				
				prot.setCoordinates(kabschFood[1]);
				wr2.write("MODEL        " + (i + 2) + "\n");
				wr2.write(prot.toString());
				wr2.write("ENDMDL\n");
			}
			wr2.close();
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}
