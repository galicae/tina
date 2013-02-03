package test;


import java.io.BufferedWriter;
import java.io.IOException;
import java.io.OutputStreamWriter;
import java.util.ArrayList;
import java.util.HashMap;

import util.JMolView;

import bioinfo.Sequence;
import bioinfo.alignment.SequenceAlignment;
import bioinfo.alignment.gotoh.FreeshiftSequenceGotoh;
import bioinfo.alignment.gotoh.GlobalSequenceGotoh;
import bioinfo.alignment.kerbsch.FreeshiftMusterLite;
import bioinfo.alignment.kerbsch.GlobalMusterLite;
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
		HashMap<String,char[]> seqlib = SeqLibrary.read("../GoBi_old/referenz/domains.seqlib");
		HashMap<String,char[]> sec_seqlib = SeqLibrary.read("secstruct.seqlib");
		ArrayList<String[]> pairs = PairReader.parse("vorescore.pairs");
		InitClass init = new InitClass();
		double[][] substMatrix = QuasarMatrix.DAYHOFF_MATRIX;
		double[][] polMatrix = init.calcGotohInputMatrix(init.calcPolarityScores());
		double[][] hbMatrix = init.calcGotohInputMatrix(init.calcHydropathyScores());
		GlobalSequenceGotoh substgotoh = new GlobalSequenceGotoh(-12.0,-1.0,substMatrix);
		FreeshiftSequenceGotoh polgotoh = new FreeshiftSequenceGotoh(-10,-4,polMatrix);
		GlobalMusterLite mustertest = new GlobalMusterLite(-20.0,-5.0,hbMatrix,polMatrix,SecStructScores.matrix,substMatrix);
			
		Sequence seq1;
		Sequence seq2;
		
		seq1 = new Sequence("1j2xA00","GPLDVQVTEDAVRRYLTRKPMTTKDLLKKFQTKKTGLSSEQTVNVLAQILKRLNPERKMINDKMHFSLK");
		seq2 = new Sequence("1p7iD00","EKRPRTAFSSEQLARLKREFNENRYLTERRRQQLSSELGLNEAQIKIWFQNARAKI");
		SequenceAlignment align1 = mustertest.align(seq1, seq2);
		SequenceAlignment align2 = substgotoh.align(seq1, seq2);
		System.out.println(align1.toStringVerbose());
		System.out.println(align2.toStringVerbose());
		
//		BufferedWriter out = new BufferedWriter(new OutputStreamWriter(System.out));
//		mustertest.streamMatricesAsTxt(out);
		
		PDBFileReader pdbreader = new PDBFileReader("../GoBi_old/STRUCTURES");
		TMMain superpos = new TMMain();
		PDBEntry p = pdbreader.readFromFolderById("1j2xA00");
		PDBEntry q = pdbreader.readFromFolderById("1p7iD00");
		Transformation tr1 = superpos.calculateTransformation(align1, p, q);
		Transformation tr2 = superpos.calculateTransformation(align2, p, q);
		System.out.println("TMScore für Muster: " + tr1.getTmscore());
		System.out.println("TMScore für Sequence: " + tr2.getTmscore());
		System.out.println("RMSD für Muster: " + tr1.getRmsd());
		System.out.println("RMSD für Sequence: " + tr2.getRmsd());
		
	
		PDBReduce reduce = new PDBReduce();
		PDBEntry[] result1 = PDBReduce.reducePDB(align1, p, tr1.transform(q));
		JMolView view = new JMolView();
		view.addPDBEntry(p, "red");
		view.addPDBEntry(q, "green");
		
		view.addPDBEntry(tr1.transform(p), "yellow");
		view.addPDBEntry(tr1.transform(q), "cyan");		
		
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
