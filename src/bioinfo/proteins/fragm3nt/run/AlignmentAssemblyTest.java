package bioinfo.proteins.fragm3nt.run;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.ByteArrayOutputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.util.LinkedList;

import org.apache.commons.exec.CommandLine;
import org.apache.commons.exec.DefaultExecutor;
import org.apache.commons.exec.PumpStreamHandler;

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
 * this class is a single entry test for the AlignmentAssembler
 * 
 * @author galicae
 * 
 */
public class AlignmentAssemblyTest {
	static int fragLength = 8;
	static int extent = 4;
	static String desktop = "/home/galicae/Desktop/STRUCTURES/";

	public static void main(String[] args) throws Exception {
		String[] bla = {"1kv3B02", "1qrkA02"};
		LinkedList<String> ids = new LinkedList<String>();
		for(int i = 0; i < bla.length; i++) {
			ids.add(bla[i]);
		}
		
		System.out.println(ids.get(0) + "\t");
		doMagic(ids);
	}

	public static void doMagic(LinkedList<String> id) throws Exception {
		// find all sequences in the id list and load them
		LinkedList<Sequence> seqs = loadSequences(id);

		// find all corresponding structures (ProteinFragment)
		LinkedList<PDBEntry> structures = loadPDBs(id);

		// now do predicting magic
		
		AlignmentAssembler ass = new AlignmentAssembler(fragLength);
		ProteinFragment pred = ass.predStrucFromAl(seqs, extent, desktop);
		PDBEntry prediction = pred.toPDB();
		double identity = pred.getClusterIndex() / 1000.0;
		System.out.println(identity + "\t");

		// now print everything in the right order: first prediction, then
		// template, then all other stuff
		BufferedWriter wr2 = new BufferedWriter(new FileWriter("fastaTest.pdb"));
		wr2.write("MODEL        1\n");
		wr2.write(pred.toString(ass.getFragments(), extent));
		wr2.write("ENDMDL\n");

		// now align every sequence with the template, kabsch structures and
		// print result
		for (int i = 0; i < 1; i++) {
			FreeshiftSequenceGotoh got = new FreeshiftSequenceGotoh(-13, -3,
					QuasarMatrix.DAYHOFF_MATRIX);
			SequenceAlignment alignment = got.align(seqs.get(0), seqs.get(i));
			PDBEntry temp = structures.get(i);
			System.out.println(temp.getID());

			double[][][] kabschFood = PDBReduce.reduce(alignment, prediction,
					temp);
			double[][][] kabschFood1 = new double[2][20][3];
			for(int j = 0; j < 20; j++) {
				kabschFood1[0][j] = kabschFood[0][j + 145];
				kabschFood1[1][j] = kabschFood[1][j + 145];
			}
			Transformation t = Kabsch.calculateTransformation(kabschFood1);
			
			System.out.println(t.getRmsd() + "\t");
			PDBEntry superposed = t.transform(temp);
			wr2.write("MODEL        " + (i + 2) + "\n");
			wr2.write(superposed.getAtomSectionAsString());
			wr2.write("ENDMDL\n");
		}
		System.out.print("\n");
		wr2.close();
		
		BufferedWriter outWr = new BufferedWriter(new FileWriter("/home/galicae/Desktop/rosenrot.pdb"));
		outWr.write(prediction.getAtomSectionAsString());
		outWr.close();
		
		String call = ("./tools/TMalign ");
		call += ("/home/galicae/Desktop/rosenrot.pdb ");
		call += (desktop + id.getFirst() + ".pdb ");
		String out = execToString(call);
		double score = findMeTmScore(out);
		System.out.println(score);
		// System.exit(0);
	}

	public static LinkedList<PDBEntry> loadPDBs(LinkedList<String> ids) {
		LinkedList<PDBEntry> pdbs = new LinkedList<PDBEntry>();
		PDBFileReader reader = new PDBFileReader(desktop);
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
	
	public static double findMeTmScore(String tm) {
		String[] r = tm.split("TM-score= ");
		String score = r[2].substring(0, 7);
		return Double.parseDouble(score);
	}

	public static String execToString(String command) throws Exception {
		ByteArrayOutputStream outputStream = new ByteArrayOutputStream();
		CommandLine commandline = CommandLine.parse(command);
		DefaultExecutor exec = new DefaultExecutor();
		PumpStreamHandler streamHandler = new PumpStreamHandler(outputStream);
		exec.setStreamHandler(streamHandler);
		exec.execute(commandline);
		return (outputStream.toString());
	}

}
