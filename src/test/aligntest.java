package test;


import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.io.OutputStreamWriter;
import java.util.ArrayList;
import java.util.HashMap;

import util.JMolView;

import bioinfo.Sequence;
import bioinfo.alignment.SequenceAlignment;
import bioinfo.alignment.gotoh.FreeshiftSequenceGotoh;
import bioinfo.alignment.gotoh.GlobalSequenceGotoh;
import bioinfo.alignment.gotoh.LocalSequenceGotoh;
import bioinfo.alignment.kerbsch.FreeshiftMusterLite;
import bioinfo.alignment.kerbsch.GLocalSequenceGotoh;
import bioinfo.alignment.kerbsch.GlobalMusterLite;
import bioinfo.alignment.kerbsch.GLocalMusterLite;
import bioinfo.alignment.kerbsch.temp.InitClass;
import bioinfo.alignment.kerbsch.temp.PairReader;
import bioinfo.alignment.kerbsch.temp.SecStructScores;
import bioinfo.alignment.kerbsch.temp.SeqLibrary;
import bioinfo.alignment.matrices.QuasarMatrix;
import bioinfo.proteins.PDBEntry;
import bioinfo.proteins.PDBFileReader;
import bioinfo.superpos.PDBReduce;
import bioinfo.superpos.TMMain;
import bioinfo.superpos.Transformation;


public class aligntest {

	public static void main(String[] args) throws Exception {
		BufferedWriter out = new BufferedWriter(new FileWriter("sec_glocal-global.aligns"));
		HashMap<String,char[]> seqlib = SeqLibrary.read("../GoBi_old/referenz/domains.seqlib");
		HashMap<String,char[]> sec_seqlib = SeqLibrary.read("secstruct.seqlib");
		ArrayList<String[]> pairs = PairReader.parse("410List.unique_pairs");
		InitClass init = new InitClass();
		double[][] secMatrix = SecStructScores.matrix;
		double[][] substMatrix = QuasarMatrix.DAYHOFF_MATRIX;
		double[][] polMatrix = init.calcGotohInputMatrix(init.calcPolarityScores());
		double[][] hbMatrix = init.calcGotohInputMatrix(init.calcHydropathyScores());
		GlobalSequenceGotoh substgotoh = new GlobalSequenceGotoh(-12.0,-1.0,substMatrix);
		FreeshiftSequenceGotoh polgotoh = new FreeshiftSequenceGotoh(-10,-4,polMatrix);
		LocalSequenceGotoh localgotoh = new LocalSequenceGotoh(-7,-2,polMatrix);
		GLocalSequenceGotoh glocalgotoh = new GLocalSequenceGotoh(-15.0,-4.0,secMatrix);
		GLocalMusterLite mustertest = new GLocalMusterLite(-4.0,-1.0, SeqLibrary.read("secstruct.seqlib"),hbMatrix,polMatrix,SecStructScores.matrix,substMatrix,0.1,0.1,0.3,0.1);
		Sequence seq1;
		Sequence seq2;
		
		String pdb1 = "1j2xA00";
		String pdb2 = "1wq2B00";
		
//		for(String[] pair : pairs){
		seq1 = new Sequence(pdb1,seqlib.get(pdb1));
		seq2 = new Sequence(pdb2,seqlib.get(pdb2));
		
		SequenceAlignment align1 = glocalgotoh.align(seq1, seq2);
		SequenceAlignment align2 = substgotoh.align(seq1, seq2);
		
		System.out.println(align1.toStringVerbose()+"\n");
		System.out.println(align2.toStringVerbose());
		
		
		PDBFileReader pdbreader = new PDBFileReader("../GoBi_old/STRUCTURES");
		TMMain superpos = new TMMain();
		PDBEntry p = pdbreader.readFromFolderById(pdb1);
		PDBEntry q = pdbreader.readFromFolderById(pdb2);
		Transformation tr1 = superpos.calculateTransformation(align1, p, q);
		Transformation tr2 = superpos.calculateTransformation(align2, p, q);
		System.out.println("TMScore für Muster: " + tr1.getTmscore());
		System.out.println("TMScore für Sequence: " + tr2.getTmscore());
		System.out.println("RMSD für Muster: " + tr1.getRmsd());
		System.out.println("RMSD für Sequence: " + tr2.getRmsd());
	
//		out.append(tr1.getRmsd()+"\t"+tr2.getRmsd()+"\n");
		
//		out.append(align1.toStringVerbose()+"\n");
//		out.append(align2.toStringVerbose()+"\n\n");
//		}
//		out.close();
	
//		PDBReduce reduce = new PDBReduce();
//		PDBEntry[] result1 = PDBReduce.reducePDB(align1, p, tr1.transform(q));
//		JMolView view = new JMolView();
//		view.addPDBEntry(p, "red");
//		view.addPDBEntry(q, "green");
//		
//		view.addPDBEntry(tr1.transform(p), "yellow");
//		view.addPDBEntry(tr1.transform(q), "cyan");		
		
//		for(String[] pair : pairs){
//			seq1 = new Sequence(pair[0],seqlib.get(pair[0]));
//			seq2 = new Sequence(pair[1],seqlib.get(pair[1]));
//			alignment = polgotoh.align(seq1, seq2);
//			System.out.println(alignment.toStringVerbose());
//		}		
		
		
//		System.out.println(alignment.toStringVerbose());
//		System.out.println(alignment.countAlignedResidues());
//		for (int i = 0; i < alignment.countAlignedResidues(); i++) {
//			System.out.print(alignment.getAlignedResidues()[0][i] + "\t");
//		}
//		System.out.println();
//		for (int i = 0; i < alignment.countAlignedResidues(); i++) {
//			System.out.print(alignment.getAlignedResidues()[1][i] + "\t");
//		}
		
	}

}
