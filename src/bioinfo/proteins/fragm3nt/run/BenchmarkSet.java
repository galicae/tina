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
import bioinfo.alignment.gotoh.GlobalSequenceGotoh;
import bioinfo.alignment.gotoh.LocalSequenceGotoh;
import bioinfo.alignment.matrices.QuasarMatrix;

/**
 * this class contains all steps of the creation of the benchmark set. This
 * comprised of 1)finding all pairs in cathscop.inpairs that belonged to the
 * same superfamily but not the same family, 2)selecting pairs with a TMalign
 * score greater than 0.5 and 3)calculating the best TM score possible
 * (comparing freeshift, global, local alignments)
 * 
 * @author galicae
 * 
 */
public class BenchmarkSet {
	public static void main(String[] args) throws Exception {
		// find all pairs from same superfamily but not same family
		// BufferedReader r = new BufferedReader(
		// new FileReader("cathscop.inpairs"));
		// String line = "";
		// LinkedList<String[]> ids = new LinkedList<String[]>();
		// while ((line = r.readLine()) != null) {
		// ids.add(line.split(" "));
		// }
		// r.close();
		// LinkedList<String[]> result = new LinkedList<String[]>();
		// System.out.println(ids.size());
		//
		// String[] l = new String[8];
		// String temp1 = l[6];
		// String temp2 = l[7];
		// String[] p1 = new String[4];
		// String[] p2 = new String[4];
		//
		// for (int i = 0; i < ids.size(); i++) {
		// l = ids.get(i);
		// if (i == 64500)
		// System.out.println("midway there");
		// temp1 = l[6];
		// temp2 = l[7];
		// p1 = temp1.split("\\.");
		// p2 = temp2.split("\\.");
		// boolean good = true;
		// for (int j = 0; j < 3; j++) {
		// if (p1[j].equals(p2[j]))
		// continue;
		// else {
		// good = false;
		// break;
		// }
		// }
		// if (good) {
		// if (!p1[3].equals(p2[3]))
		// result.add(l);
		// }
		// }
		// System.out.println(result.size()
		// + " pairs in the same superfamily but not same family");
		// BufferedWriter w = new BufferedWriter(new FileWriter("spnf"));
		// for (String[] s : result) {
		// w.write(s[0] + " " + s[1] + "\n");
		// }
		// w.close();

		// find all pairs with a tm score larger than 0.5
		// BufferedReader r = new BufferedReader(new FileReader("spnf"));
		// BufferedWriter w = new BufferedWriter(new FileWriter("spnfTM"));
		// String line = "";
		// String alignOut = "";
		// StringBuilder call = new StringBuilder();
		// String desktop = "/home/galicae/Desktop/STRUCTURES/";
		// String[] pair = new String[2];
		// double tempTM = 0;
		//
		// while ((line = r.readLine()) != null) {
		// pair = line.split(" ");
		// System.out.println(line);
		// call.setLength(0);
		// call.append("./lib/TMalign ");
		// call.append(desktop + pair[0] + ".pdb ");
		// call.append(desktop + pair[1] + ".pdb -a");
		// alignOut = execToString(call.toString());
		// tempTM = findMeTmScore(alignOut);
		// if (tempTM >= 0.5) {
		// w.write(line + " " + tempTM + "\n");
		// System.err.println("got a pair");
		// }
		// }
		// r.close();
		// w.close();

		// calculate gotoh TM score for said alignments
		BufferedReader r = new BufferedReader(new FileReader("spnfTM"));

		String line = "";
		FreeshiftSequenceGotoh f = new FreeshiftSequenceGotoh(-9, -3,
				QuasarMatrix.DAYHOFF_MATRIX);
		GlobalSequenceGotoh g = new GlobalSequenceGotoh(-12, -1,
				QuasarMatrix.DAYHOFF_MATRIX);
		LocalSequenceGotoh l = new LocalSequenceGotoh(-12, -1,
				QuasarMatrix.DAYHOFF_MATRIX);

		double fScore = 0;
		double gScore = 0;
		double lScore = 0;

		String fcall = "";
		String gcall = "";
		String lcall = "";

		String fOut = "";
		String gOut = "";
		String lOut = "";

		LinkedList<Double> gotohPerformance = new LinkedList<Double>();

		String desktop = "/home/galicae/Desktop/STRUCTURES/";

		while ((line = r.readLine()) != null) {
			// calculate all possible alignments to find the one with the best
			// TMalign score

			String[] arr = line.split(" ");
			Sequence seq1 = retrieveSeq(arr[0]);
			Sequence seq2 = retrieveSeq(arr[1]);

			SequenceAlignment fAli = f.align(seq1, seq2);
			SequenceAlignment gAli = g.align(seq1, seq2);
			SequenceAlignment lAli = l.align(seq1, seq2);

			BufferedWriter lw = new BufferedWriter(new FileWriter(
					"./fastaFiles/loc" + seq1.getID() + "_" + seq2.getID()));
			BufferedWriter fw = new BufferedWriter(new FileWriter(
					"./fastaFiles/fre" + seq1.getID() + "_" + seq2.getID()));
			BufferedWriter gw = new BufferedWriter(new FileWriter(
					"./fastaFiles/glo" + seq1.getID() + "_" + seq2.getID()));
			lw.write(toFastaFormat(lAli));
			fw.write(toFastaFormat(fAli));
			gw.write(toFastaFormat(gAli));
			lw.close();
			fw.close();
			gw.close();

			// try finding the best score
			lcall = ("./lib/TMalign ");
			lcall += (desktop + seq1.getID() + ".pdb ");
			lcall += (desktop + seq2.getID() + ".pdb ");
			lcall += "-I ./fastaFiles/loc" + seq1.getID() + "_" + seq2.getID();
			lOut = execToString(lcall);
			lScore = findMeTmScore(lOut);

			gcall = ("./lib/TMalign ");
			gcall += (desktop + seq1.getID() + ".pdb ");
			gcall += (desktop + seq2.getID() + ".pdb ");
			gcall += "-I ./fastaFiles/glo" + seq1.getID() + "_" + seq2.getID();
			gOut = execToString(gcall);
			gScore = findMeTmScore(gOut);

			fcall = ("./lib/TMalign ");
			fcall += (desktop + seq1.getID() + ".pdb ");
			fcall += (desktop + seq2.getID() + ".pdb ");
			fcall += "-I ./fastaFiles/fre" + seq1.getID() + "_" + seq2.getID();
			fOut = execToString(fcall);
			fScore = findMeTmScore(fOut);

			gotohPerformance.add(Math.max(Math.max(fScore, lScore), gScore));
		}

		BufferedWriter tmw = new BufferedWriter(new FileWriter("gotohTMScore"));
		for (double d : gotohPerformance) {
			tmw.write(d + "\n");
		}
		tmw.close();

		r.close();
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

	public static Sequence retrieveSeq(String id) throws Exception {
		BufferedReader r = new BufferedReader(new FileReader("domains.seqlib"));
		String line = "";
		while ((line = r.readLine()) != null) {
			if (line.startsWith(id)) {
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
