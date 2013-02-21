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
import bioinfo.alignment.kerbsch.FBGotoh;
import bioinfo.alignment.kerbsch.FreeshiftMusterLite;
import bioinfo.alignment.kerbsch.GLocalSequenceGotoh;
import bioinfo.alignment.kerbsch.GlobalMusterLite;
import bioinfo.alignment.kerbsch.GLocalMusterLite;
import bioinfo.alignment.kerbsch.HSPAlignment;
import bioinfo.alignment.kerbsch.Kerbsch;
import bioinfo.alignment.kerbsch.temp.InitClass;
import bioinfo.alignment.kerbsch.temp.PairReader;
import bioinfo.alignment.kerbsch.temp.SecStructScores;
import bioinfo.alignment.kerbsch.temp.SeqLibrary;
import bioinfo.alignment.matrices.QuasarMatrix;
import bioinfo.proteins.PDBEntry;
import bioinfo.proteins.PDBFileReader;
import bioinfo.superpos.Kabsch;
import bioinfo.superpos.PDBReduce;
import bioinfo.superpos.TMMain;
import bioinfo.superpos.Transformation;


public class aligntest {

	public static void main(String[] args) throws Exception {
		BufferedWriter out = null;
		out = new BufferedWriter(new FileWriter("pol_local.scores"));
		HashMap<String,char[]> seqlib = SeqLibrary.read("../full_domains.seqlib");
		HashMap<String,char[]> sec_seqlib = SeqLibrary.read("../secstruct.seqlib");
		ArrayList<String[]> pairs = PairReader.parse("../410List.unique_pairs");
		InitClass init = new InitClass();
//		double[][] secMatrix = SecStructScores.matrix;
//		double[][] substMatrix = QuasarMatrix.DAYHOFF_MATRIX;
//		double[][] polMatrix = init.calcGotohInputMatrix(init.calcPolarityScores());
//		double[][] hbMatrix = init.calcGotohInputMatrix(init.calcHydropathyScores());
		Kerbsch kerbsch = new Kerbsch(-20,-4,QuasarMatrix.DAYHOFF_MATRIX,init.calcGotohInputMatrix(init.calcHydropathyScores()),init.calcGotohInputMatrix(init.calcPolarityScores()),SecStructScores.matrix,sec_seqlib);

//		GlobalSequenceGotoh substgotoh = new GlobalSequenceGotoh(-12.0,-1.0,substMatrix);
//		FreeshiftSequenceGotoh polgotoh = new FreeshiftSequenceGotoh(-10,-4,polMatrix);
		LocalSequenceGotoh localgotoh = new LocalSequenceGotoh(-12.0,-1.0,QuasarMatrix.DAYHOFF_MATRIX);
//		GLocalSequenceGotoh glocalgotoh = new GLocalSequenceGotoh(-15.0,-4.0,secMatrix);
//		GlobalMusterLite mustertest = new GlobalMusterLite(-6.0,-1.0, SeqLibrary.read("../secstruct.seqlib"),init.calcGotohInputMatrix(init.calcHydropathyScores()),init.calcGotohInputMatrix(init.calcPolarityScores()),SecStructScores.matrix,substMatrix,0.0,0.1,0.5,0.4);
//		HSPAlignment hsptest = new HSPAlignment(-12,-1,2,2,2,QuasarMatrix.DAYHOFF_MATRIX);
		FBGotoh fbgotoh = new FBGotoh(-12.0,-1.0,0.38,6.8,61735,QuasarMatrix.DAYHOFF_MATRIX);
		Sequence seq1;
		Sequence seq2;
		BufferedWriter writer = new BufferedWriter(new OutputStreamWriter(System.out));
		String pdb1 = "1tuaA02";
		String pdb2 = "1tuaA01";
		
//		for(String[] pair : pairs){	
//			pdb1 = pair[0];
//			pdb2 = pair[1];
//			seq1 = new Sequence(pdb1,seqlib.get(pdb1));
//			seq2 = new Sequence(pdb2,seqlib.get(pdb2));
//			kerbsch.align(seq1, seq2);
//		}
		
		seq1 = new Sequence(pdb1,seqlib.get(pdb1));
		seq2 = new Sequence(pdb2,seqlib.get(pdb2));
		SequenceAlignment kalign = kerbsch.align(seq1, seq2);
		System.out.println(kalign.toStringVerbose());

		System.out.println(localgotoh.align(seq1, seq2).toStringVerbose());
//		kerbsch.streamMatricesAsTxt(writer);
		
//		SequenceAlignment align2 = mustertest.align(seq1, seq2);
//		System.out.println(align1.toStringVerbose()+"\n");
//		System.out.println(align2.toStringVerbose()+"\n");
		
//		PDBFileReader pdbreader = new PDBFileReader("../GoBi_old/STRUCTURES");
//		PDBEntry p = pdbreader.readFromFolderById(pdb1);
//		PDBEntry q = pdbreader.readFromFolderById(pdb2);
//		Transformation tr1 = Kabsch.calculateTransformation(PDBReduce.reduce(align1, p, q));
//		Transformation tr2 = Kabsch.calculateTransformation(PDBReduce.reduce(align2, p, q));
//
//		System.out.println("TMScore for Kerbsch: " + tr1.getTmscore());
//		System.out.println("TMScore for muster: " + tr2.getTmscore());
//		System.out.println("RMSD for Kerbsch: " + tr1.getRmsd());
//		System.out.println("RMSD for muster: " + tr2.getRmsd());
		

//		TMMain superpos1 = new TMMain();
//		TMMain superpos12 = new TMMain();

//		Transformation tr1 = superpos1.calculateTransformation(align1, p, q);
//		Transformation tr2 = superpos12.calculateTransformation(align2, p, q);

//	
//		out.append(tr1.getRmsd()+"\t"+tr2.getRmsd()+"\n");
		
//		out.append(align1.toStringVerbose()+"\n");
//		out.append(align2.toStringVerbose()+"\n\n");
//		}
//		out.close();
//	
//
//		JMolView view = new JMolView();
//		
//		view.addPDBEntry(p, "yellow");
//		view.addPDBEntry(tr2.transform(q), "cyan");		
		
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
