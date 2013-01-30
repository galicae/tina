package test;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.InputStreamReader;

import bioinfo.alignment.SequenceAlignment;
import bioinfo.alignment.SequenceAlignmentFileReader;
import bioinfo.proteins.PDBEntry;
import bioinfo.proteins.PDBFileReader;
import bioinfo.superpos.TMCollective;
import bioinfo.superpos.TMMain;
import bioinfo.superpos.TMOriginal;
import bioinfo.superpos.Transformation;

public class TMPipeline {
	public static void main(String[] args) throws Exception {
		String protDir = "/home/p/papadopoulos/Desktop/STRUCTURES/";
		SequenceAlignmentFileReader aliReader = new SequenceAlignmentFileReader(
				"/home/p/papadopoulos/Dokumente/assignment1/out/freeshift.alignments");
		aliReader.initSequentialRead();

		PDBEntry[] tmOr = new PDBEntry[2];

		TMCollective tmc = new TMCollective();
		PDBFileReader pdbReader = new PDBFileReader();
		SequenceAlignment curAli = aliReader.nextAlignment();
		curAli = aliReader.nextAlignment();
		curAli = aliReader.nextAlignment();

		PDBEntry pdb1 = pdbReader.readPDBFromFile(protDir
				+ curAli.getComponent(0).getID() + ".pdb");
		PDBEntry pdb2 = pdbReader.readPDBFromFile(protDir
				+ curAli.getComponent(1).getID() + ".pdb");

		TMMain main = new TMMain();
		System.out.println(curAli.toStringVerbose() + "\n");
		System.err.println("OBJECT ORIENTED#########################");
		Transformation tr1 = main.calculateTransformation(curAli, pdb1, pdb2);
		System.out.println(tr1.getTmscore());

		tmOr = tmc.createTMInput(curAli, pdb1, pdb2);

		String p2 = "TM" + curAli.getComponent(1).getID() + ".pdb";
		String p1 = "TM" + curAli.getComponent(0).getID() + ".pdb";
		try {
			BufferedWriter wr1 = new BufferedWriter(new FileWriter(p1));
			wr1.write(tmOr[0].getAtomSectionAsString());
			wr1.close();
			BufferedWriter wr2 = new BufferedWriter(new FileWriter(p2));
			wr2.write(tmOr[1].getAtomSectionAsString());
			wr2.close();
			// Runtime r = Runtime.getRuntime();
			// Process p = r.exec("./lib/TMscore " + p1 + " " + p2);
			// + " -l "
			// + curAli.getComponent(1).length());
			// BufferedReader input = new BufferedReader(new InputStreamReader(
			// p.getInputStream()));
			//
			// String l = input.readLine();
			// while (l != null) {
			// // if(l.startsWith("TM-score    ="))
			// System.out.println(l);
			// l = input.readLine();
			// }
			//
			// p.destroy();
		} catch (Exception e) {
			e.printStackTrace();
		}
		System.err.println("REFERENCE#################################");
		TMOriginal.calculateTmScore(p1, p2, curAli.getComponent(1).length());

		// while (curAli != null) {
		// pdb1 = pdbReader.readPDBFromFile(protDir
		// + curAli.getComponent(0).getID() + ".pdb");
		// pdb2 = pdbReader.readPDBFromFile(protDir
		// + curAli.getComponent(1).getID() + ".pdb");
		// tr1 = main.calculateTransformation(curAli, pdb1, pdb2);
		// System.out.println(tr1.getTmscore());
		// curAli = aliReader.nextAlignment();
		// }
	}
}
