package bioinfo.proteins.fragm3nt.run;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.util.HashMap;
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
import bioinfo.superpos.TMMain;
import bioinfo.superpos.Transformation;

public class AlignmentAssemblyRun {
	static int fragLength = 8;
	static String pdbDirectory = "/home/galicae/Desktop/STRUCTURES/";
	static BufferedWriter resultWriter;

	public static void main(String[] args) throws Exception {
		resultWriter  = new BufferedWriter(new FileWriter("bucket05Results"));
		// bucket file direction
		String bucket = "bucket05";// args[0];
		// first find how many different sequences there are in the bucket
		// so first read bucket file
		// idea: use every ID as key and map to it a list with all the sequence
		// IDs that are aligned against the first ID
		HashMap<String, LinkedList<String>> map = new HashMap<String, LinkedList<String>>();
		BufferedReader r = new BufferedReader(new FileReader(bucket));

		String line = "";
		String[] lineArr = new String[2];
		while ((line = r.readLine()) != null) {
			lineArr = line.split(" ");
			if (!map.containsKey(lineArr[0])) {
				map.put(lineArr[0], new LinkedList<String>());
				map.get(lineArr[0]).add(lineArr[0]);
				map.get(lineArr[0]).add(lineArr[1]);
			} else {
				map.get(lineArr[0]).add(lineArr[1]);
			}
			if (!map.containsKey(lineArr[1])) {
				map.put(lineArr[1], new LinkedList<String>());
				map.get(lineArr[1]).add(lineArr[1]);
				map.get(lineArr[1]).add(lineArr[0]);
			} else {
				map.get(lineArr[1]).add(lineArr[0]);
			}
		}
		r.close();

		// now for every entry do the stuff from AlignmentAssemblyTest
		for (String s : map.keySet()) {
			try {
				resultWriter.write(s + "\t");
				System.out.println(s + "\t");
				doMagic(map.get(s));
			} catch (Exception e) {
				e.printStackTrace();
				continue;
			}
		}
		resultWriter.close();
	}

	public static void doMagic(LinkedList<String> id) throws Exception {
		// find all sequences in the id list and load them
		LinkedList<Sequence> seqs = loadSequences(id);

		// find all corresponding structures (ProteinFragment)
		LinkedList<PDBEntry> structures = loadPDBs(id);

		// now do predicting magic
		int extent = 5;
		AlignmentAssembler ass = new AlignmentAssembler(fragLength);
		ProteinFragment pred = ass.predStrucFromAl(seqs, extent, pdbDirectory);
		PDBEntry prediction = pred.toPDB();
		double identity = pred.getClusterIndex() / 1000.0;
		resultWriter.write(identity + "\t");

		// now print everything in the right order: first prediction, then
		// template, then all other stuff
		BufferedWriter wr2 = new BufferedWriter(new FileWriter("./bucket5/" + id.get(0) + ".pdb"));
		wr2.write("MODEL        1\n");
		wr2.write(pred.toString(ass.getFragments(), extent));
		wr2.write("ENDMDL\n");
		wr2.write("MODEL        2\n");
		wr2.write(structures.get(0).getAtomSectionAsString());
		wr2.write("ENDMDL\n");
		// now align every sequence with the template, kabsch structures and print result
		for(int i = 1; i < seqs.size(); i++) {
			FreeshiftSequenceGotoh got = new FreeshiftSequenceGotoh(-13, -3, QuasarMatrix.DAYHOFF_MATRIX);
			SequenceAlignment alignment = got.align(seqs.get(0), seqs.get(i));
			PDBEntry temp = structures.get(i);
			
			double[][][] kabschFood = PDBReduce.reduce(alignment, structures.get(0), temp);
			Transformation t = Kabsch.calculateTransformation(kabschFood);
			resultWriter.write(t.getRmsd() + "\t");
			PDBEntry superposed = t.transform(temp);
			wr2.write("MODEL        " + (i + 2) + "\n");
			wr2.write(superposed.getAtomSectionAsString());
			wr2.write("ENDMDL\n");
		}
		resultWriter.write("\n");
		wr2.close();
//		System.exit(0);
	}

	public static LinkedList<PDBEntry> loadPDBs(LinkedList<String> ids) {
		LinkedList<PDBEntry> pdbs = new LinkedList<PDBEntry>();
		PDBFileReader reader = new PDBFileReader(pdbDirectory);
		for (int i = 0; i < ids.size(); i++) {
			pdbs.add(reader.readFromFolderById(ids.get(i)));
		}
		return pdbs;
	}

	public static LinkedList<Sequence> loadSequences(LinkedList<String> ids) {
		LinkedList<Sequence> result = new LinkedList<Sequence>();
		try {
			BufferedReader r = new BufferedReader(new FileReader(
					"domains.seqlib"));
			String line = "";
			int count = 0;
			String stringSeq = "";
			while (count < ids.size() && (line = r.readLine()) != null) {
				for (int i = 0; i < ids.size(); i++) {
					if (line.startsWith(ids.get(i))) {
						stringSeq = line.split(":")[1];
						result.add(new Sequence(ids.get(i), stringSeq));
						count++;
						continue;
					}
				}
			}
			r.close();
		} catch (Exception e) {
			e.printStackTrace();
		}
		return result;
	}

}
