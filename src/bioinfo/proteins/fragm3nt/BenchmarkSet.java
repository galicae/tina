package bioinfo.proteins.fragm3nt;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.ByteArrayOutputStream;
import java.io.FileReader;
import java.io.FileWriter;

import org.apache.commons.exec.CommandLine;
import org.apache.commons.exec.DefaultExecutor;
import org.apache.commons.exec.PumpStreamHandler;

public class BenchmarkSet {
	public static void main(String[] args) throws Exception {
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
		BufferedReader r = new BufferedReader(new FileReader("spnf"));
		BufferedWriter w = new BufferedWriter(new FileWriter("spnfTM"));
		String line = "";
		String alignOut = "";
		StringBuilder call = new StringBuilder();
		String desktop = "/home/galicae/Desktop/STRUCTURES/";
		String[] pair = new String[2];
		double tempTM = 0;

		while ((line = r.readLine()) != null) {
			pair = line.split(" ");
			System.out.println(line);
			call.setLength(0);
			call.append("./lib/TMalign ");
			call.append(desktop + pair[0] + ".pdb ");
			call.append(desktop + pair[1] + ".pdb -a");
			alignOut = execToString(call.toString());
			tempTM = findMeTmScore(alignOut);
			if (tempTM >= 0.5) {
				w.write(line + " " + tempTM + "\n");
				System.err.println("got a pair");
			}
		}
		r.close();
		w.close();
	}

	public static double findMeTmScore(String tm) {
		String[] r = tm.split("TM-score= ");
		String score = r[3].substring(0, 7);
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
