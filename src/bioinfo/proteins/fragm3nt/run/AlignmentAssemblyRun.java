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
		
		// now align every sequence with the template, kabsch structures and print result
		for(int i = 0; i < seqs.size(); i++) {
			FreeshiftSequenceGotoh got = new FreeshiftSequenceGotoh(-13, -3, QuasarMatrix.DAYHOFF_MATRIX);
			SequenceAlignment alignment = got.align(seqs.get(0), seqs.get(i));
			PDBEntry temp = structures.get(i);
			
			double[][][] kabschFood = PDBReduce.reduce(alignment, prediction, temp);
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
//		for (int i = 0; i < ids.size(); i++) {
//			ProteinFragment temp = new ProteinFragment(ids.get(i),
//					PDBReduce.reduceSinglePDB(pdbs.get(i)), fragLength);
//			result.add(temp);
//		}
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

	// public static void doMagic(LinkedList<String> seqIdList) {
	// String[] seqIds = new String[seqIdList.size()];
	// for (int i = 0; i < seqIds.length; i++) {
	// seqIds[i] = seqIdList.get(i);
	// }
	// String[] stringSeq = new String[seqIds.length];
	//
	// // read sequences (strings)
	// try {
	// BufferedReader r = new BufferedReader(new FileReader(
	// "domains.seqlib"));
	// String line = "";
	// int count = 0;
	// while (count < seqIds.length && (line = r.readLine()) != null) {
	// for (int i = 0; i < seqIds.length; i++) {
	// if (line.startsWith(seqIds[i])) {
	// stringSeq[i] = line.split(":")[1];
	// count++;
	// continue;
	// }
	// }
	// }
	// r.close();
	// } catch (Exception e) {
	// e.printStackTrace();
	// }
	//
	// LinkedList<Sequence> seqs = new LinkedList<Sequence>();
	// for (int i = 0; i < stringSeq.length; i++) {
	// seqs.add(new Sequence(seqIds[i], stringSeq[i]));
	// }
	//
	// int extent = 5;
	// AlignmentAssembler ass = new AlignmentAssembler(8);
	// ProteinFragment result = ass.predStrucFromAl(seqs, extent,
	// "/home/galicae/Desktop/STRUCTURES/");
	// PDBFileReader read = new PDBFileReader(
	// "/home/galicae/Desktop/STRUCTURES/");
	// PDBEntry pdb = read.readFromFolderById(seqIds[0]);
	// String query = "";
	// for (int i = 0; i < pdb.length(); i++) {
	// query += pdb.getAminoAcid(i).getName();
	// }
	// double[][][] kabschFood = new double[2][][];
	//
	// try {
	// BufferedWriter wr2 = new BufferedWriter(new FileWriter("./bucket5/"
	// + seqIds[0] + ".pdb"));
	// kabschFood = new double[2][result.getAllResidues().length][3];
	// kabschFood[0] = result.getAllResidues();
	// kabschFood[1] = PDBReduce.reduceSinglePDB(pdb);
	// Transformation t = Kabsch.calculateTransformation(kabschFood);
	// kabschFood[1] = t.transform(kabschFood[1]);
	// resultWriter.writeln(t.getRmsd());
	//wr2.write("MODEL        " + (i + 2) + "\n");
	// ProteinFragment prot = new ProteinFragment("real", "D",
	// new double[0][0], 8);
	// prot.setSequence(query);
	// prot.append(PDBReduce.reduceSinglePDB(pdb), "");
	// pdb = t.transform(pdb);
	// prot.setCoordinates(kabschFood[1]);
	//
	// wr2.write("MODEL        1\n");
	// wr2.write(result.toString(ass.getFragments(), extent));
	// wr2.write("ENDMDL\n");
	// wr2.write("MODEL        2\n");
	// wr2.write(prot.toString());
	// wr2.write("ENDMDL\n");
	//
	// // kabschFood[0] = prot.getAllResidues();
	// Sequence seq1 = new Sequence(pdb.getID(), pdb.getSequence());
	// for (int i = 1; i < seqs.size(); i++) {
	// PDBEntry pdb1 = read.readFromFolderById(seqs.get(i).getID());
	// FreeshiftSequenceGotoh got = new FreeshiftSequenceGotoh(-15,
	// -3, QuasarMatrix.DAYHOFF_MATRIX);
	// Sequence seq2 = new Sequence(pdb1.getID(), pdb1.getSequence());
	// SequenceAlignment ali = got.align(seq1, seq2);
	// // resultWriter.writeln(ali.toStringVerbose());
	//
	// prot = new ProteinFragment(pdb1.getID(), "D", new double[0][0],
	// 8);
	// prot.setSequence(pdb1.getSequence());
	// prot.append(PDBReduce.reduceSinglePDB(pdb1), "");
	//
	// kabschFood[1] = prot.getAllResidues();
	// kabschFood = PDBReduce.reduce(ali, pdb, pdb1);
	// t = Kabsch.calculateTransformation(kabschFood);
	// kabschFood[1] = t.transform(kabschFood[1]);
	//
	// prot.setCoordinates(kabschFood[1]);
	// wr2.write("MODEL        " + (i + 2) + "\n");
	// wr2.write(prot.toString());
	// wr2.write("ENDMDL\n");
	// }
	// wr2.close();
	// } catch (Exception e) {
	// e.printStackTrace();
	// }
	// }

}
