package bioinfo.proteins.fragm3nt.run;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.ByteArrayOutputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.util.HashMap;
import java.util.LinkedList;

import org.apache.commons.exec.CommandLine;
import org.apache.commons.exec.DefaultExecutor;
import org.apache.commons.exec.PumpStreamHandler;

import bioinfo.Sequence;
import bioinfo.alignment.SequenceAlignment;
import bioinfo.alignment.gotoh.FreeshiftSequenceGotoh;
import bioinfo.alignment.gotoh.GlobalSequenceGotoh;
import bioinfo.alignment.gotoh.LocalSequenceGotoh;
import bioinfo.alignment.matrices.QuasarMatrix;
import bioinfo.proteins.PDBEntry;
import bioinfo.proteins.PDBFileReader;
import bioinfo.proteins.fragm3nt.ProteinFragment;
import bioinfo.proteins.fragm3nt.assembler.AlignmentAssembler;
import bioinfo.superpos.TMMain;
import bioinfo.superpos.Transformation;

/**
 * purpose of this class is to evaluate the TM score of the models produced by
 * the AlignmentAssembler class on the bucket05 data
 * 
 * @author galicae
 * 
 */
public class AliAssVsGotoh {
	static int fragLength = 9;
	static int extent = 4;
	static String desktop = "/home/galicae/Desktop/STRUCTURES/";
	static BufferedWriter resultWriter;

	public static void main(String[] args) throws Exception {
		BufferedWriter resultWriter  = new BufferedWriter(new FileWriter("gotohVSass"));
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
			System.out.println(s);
			try {
				// calculate Gotoh score
				resultWriter.write(s + "\t");
				
				String seq1ID = map.get(s).get(0);
				String seq2ID = map.get(s).get(1);
				double gotScore = gotohPart(seq1ID, seq2ID);
				System.out.print(gotScore + "\t");
				resultWriter.write(gotScore + "\t");
				// calculate Fragm3nt score
				
				double assScore = doMagic(map.get(s));
				System.out.println(assScore);
				resultWriter.write(assScore + "\n");
			} catch (Exception e) {
				e.printStackTrace();
//				System.exit(0);
				continue;
			}
		}
		resultWriter.close();
	}
	
public static double gotohPart(String seq1ID, String seq2ID) throws Exception {
		
		FreeshiftSequenceGotoh f = new FreeshiftSequenceGotoh(-9, -3, QuasarMatrix.DAYHOFF_MATRIX);
		GlobalSequenceGotoh g = new GlobalSequenceGotoh(-12, -1, QuasarMatrix.DAYHOFF_MATRIX);
		LocalSequenceGotoh l = new LocalSequenceGotoh(-12, -1, QuasarMatrix.DAYHOFF_MATRIX);
		
		Sequence seq1 = retrieveSeq(seq1ID);
		Sequence seq2 = retrieveSeq(seq2ID);
		
		SequenceAlignment fAli = f.align(seq1, seq2);
		SequenceAlignment gAli = g.align(seq1, seq2);
		SequenceAlignment lAli = l.align(seq1, seq2);
		
		BufferedWriter lw = new BufferedWriter(new FileWriter("./fastaFiles/loc" + seq1.getID() + "_" + seq2.getID()));
		BufferedWriter fw = new BufferedWriter(new FileWriter("./fastaFiles/fre" + seq1.getID() + "_" + seq2.getID()));
		BufferedWriter gw = new BufferedWriter(new FileWriter("./fastaFiles/glo" + seq1.getID() + "_" + seq2.getID()));
//		lw.write(toFastaFormat(lAli));
//		fw.write(toFastaFormat(fAli));
//		gw.write(toFastaFormat(gAli));
		lw.write(lAli.toStringVerbose());
		fw.write(fAli.toStringVerbose());
		gw.write(gAli.toStringVerbose());
		lw.close();
		fw.close();
		gw.close();
		
		double fScore = 0;
		double gScore = 0;
		double lScore = 0;
		
		String fcall = "";
		String gcall = "";
		String lcall = "";
		
		String fOut = "";
		String gOut = "";
		String lOut = "";
		
//		lcall = ("./tools/TMalign ");
//		lcall += (desktop + seq1ID + ".pdb ");
//		lcall += (desktop + seq2ID + ".pdb ");
//		lcall += "-I ./fastaFiles/loc" + seq1ID + "_" + seq2ID;
//		lOut = execToString(lcall);
//		lScore = findMeTmScore(lOut);
//		
//		gcall = ("./tools/TMalign ");
//		gcall += (desktop + seq1ID + ".pdb ");
//		gcall += (desktop + seq2ID + ".pdb ");
//		gcall += "-I ./fastaFiles/glo" + seq1ID + "_" + seq2ID;
//		gOut = execToString(gcall);
//		gScore = findMeTmScore(gOut);
//		
//		fcall = ("./tools/TMalign ");
//		fcall += (desktop + seq1ID + ".pdb ");
//		fcall += (desktop + seq2ID + ".pdb ");
//		fcall += "-I ./fastaFiles/fre" + seq1ID + "_" + seq2ID;
//		fOut = execToString(fcall);
//		fScore = findMeTmScore(fOut);
		
		TMMain m = new TMMain();
		String fLoc = "./fastaFiles/fre" + seq1.getID() + "_" + seq2.getID();
		String gLoc = "./fastaFiles/glo" + seq1.getID() + "_" + seq2.getID();
		String lLoc = "./fastaFiles/loc" + seq1.getID() + "_" + seq2.getID();
		
		String pdb1 = desktop + seq1ID + ".pdb";
		String pdb2 = desktop + seq2ID + ".pdb";
		
		int nLeng = fAli.getComponent(1).length();
		
		Transformation ft = m.calculateTransformation(fLoc, pdb1, pdb2, nLeng);
		Transformation gt = m.calculateTransformation(gLoc, pdb1, pdb2, nLeng);
		Transformation lt = m.calculateTransformation(lLoc, pdb1, pdb2, nLeng);
		
		fScore = ft.getTmscore();
		gScore = gt.getTmscore();
		lScore = lt.getTmscore();
		
		return Math.max(fScore, Math.max(gScore, lScore));
	}

	public static double doMagic(LinkedList<String> id) throws Exception {
		// find all sequences in the id list and load them
		LinkedList<Sequence> seqs = loadSequences(id);

		// find all corresponding structures (ProteinFragment)
		LinkedList<PDBEntry> structures = loadPDBs(id);

		// now do predicting magic
		AlignmentAssembler ass = new AlignmentAssembler(fragLength);
		ProteinFragment pred = ass.predStrucFromAl(seqs, extent, desktop);
		PDBEntry prediction = pred.toPDB();
		
//		TMMain main = new TMMain();
//		Sequence s = seqs.getFirst();
//		char[] row = s.getSequence();
//		SequenceAlignment alignment = new SequenceAlignment(s, s, row, row, 5);
//		Transformation t = main.calculateTransformation(alignment, structures.get(0), prediction);
//		return t.getTmscore();
		
		// write pdb to file
		BufferedWriter write5r = new BufferedWriter(new FileWriter("predict.pdb"));
		write5r.write(prediction.getAtomSectionAsString());
		write5r.close();
		String call = ("./tools/TMalign ");
		call += ("./predict.pdb ");
		call += (desktop + id.get(0) + ".pdb");
		String out = execToString(call);
		double assScore = findMeTmScore(out);
		
		return assScore;
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
		String score = r[1].substring(0, 7);
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
	
	public static Sequence retrieveSeq(String id) throws Exception {
		BufferedReader r = new BufferedReader(new FileReader("domains.seqlib"));
		String line = "";
		while((line = r.readLine()) != null) {
			if(line.startsWith(id)){
				String[] s = line.split(":");
				Sequence sq = new Sequence(s[0], s[1]);
				r.close();
				return sq;
			}
		}
		r.close();
		return null;
	}
	
	public static String toFastaFormat(SequenceAlignment ali) {
		StringBuilder sb = new StringBuilder();
		sb.append(">" + ali.getComponent(0).getID() + "\n");
		for (int i = 0; i < ali.getRowAsString(0).length(); i++) {
			sb.append(ali.getRow(0)[i]);
			if (i % 80 == 0 && i > 0)
				sb.append("\n");
		}
		sb.append("\n");
		sb.append(">" + ali.getComponent(1).getID() + "\n");
		for (int i = 0; i < ali.getRowAsString(1).length(); i++) {
			sb.append(ali.getRow(1)[i]);
			if (i % 80 == 0 && i > 0)
				sb.append("\n");
		}
		return sb.toString();
	}
	
	
}
