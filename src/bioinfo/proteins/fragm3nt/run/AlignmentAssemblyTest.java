package bioinfo.proteins.fragm3nt.run;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.util.LinkedList;

import bioinfo.Sequence;
import bioinfo.proteins.PDBEntry;
import bioinfo.proteins.PDBFileReader;
import bioinfo.proteins.fragm3nt.ProteinFragment;
import bioinfo.proteins.fragm3nt.assembler.AlignmentAssembler;
import bioinfo.superpos.Kabsch;
import bioinfo.superpos.PDBReduce;
import bioinfo.superpos.Transformation;

public class AlignmentAssemblyTest {

	public static void main(String[] args) {
		String[] seqIds = { "1f0yA02", "3hdhC02"};
		String[] stringSeq = new String[seqIds.length];
		
		// read sequences (strings)
		try {
			BufferedReader r = new BufferedReader(new FileReader("domains.seqlib"));
			String line = "";
			int count = 0;
			while(count < seqIds.length && (line=r.readLine()) != null) {
				for(int i = 0; i < seqIds.length; i++) {
					if(line.startsWith(seqIds[i])) {
						stringSeq[i] = line.split(":")[1];
						count++;
						continue;
					}
				}
			}
			r.close();
		}
		catch(Exception e) {
			e.printStackTrace();
		}
		
		LinkedList<Sequence> seqs = new LinkedList<Sequence>();
		for(int i = 0; i < stringSeq.length; i++) {
			seqs.add(new Sequence(seqIds[i], stringSeq[i]));
		}

		AlignmentAssembler ass = new AlignmentAssembler(8);
		ProteinFragment result = ass.predStrucFromAl(seqs, 5, "/home/galicae/Desktop/STRUCTURES/");
		PDBFileReader read = new PDBFileReader("/home/galicae/Desktop/STRUCTURES/");
		PDBEntry pdb = read.readFromFolderById("1j2xA00");
		String query = "";
		for(int i = 0; i < pdb.length(); i++) {
			query += pdb.getAminoAcid(i).getName();
		}
		double[][][] kabschFood = new double[2][][];
		
		try {
			BufferedWriter wr2 = new BufferedWriter(new FileWriter("fastaTest.pdb"));
			kabschFood = new double[2][result.getAllResidues().length][3];
			kabschFood[0] = result.getAllResidues();
			kabschFood[1] = PDBReduce.reduceSinglePDB(pdb);
			Transformation t = Kabsch
					.calculateTransformation(kabschFood);
			kabschFood[1] = t.transform(kabschFood[1]);
			System.out.println(t.getRmsd());

			ProteinFragment prot = new ProteinFragment("real", "D", new double[0][0], 8);
			prot.setSequence(query);
			prot.append(PDBReduce.reduceSinglePDB(pdb), "");
			prot.setCoordinates(kabschFood[1]);

			wr2.write("MODEL        1\n");
			wr2.write(result.toString());
			wr2.write("ENDMDL\n");
			wr2.write("MODEL        2\n");
			wr2.write(prot.toString());
			wr2.write("ENDMDL\n");
			wr2.close();
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}
