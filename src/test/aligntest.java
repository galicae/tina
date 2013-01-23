package test;


import java.util.ArrayList;
import java.util.HashMap;

import bioinfo.Sequence;
import bioinfo.alignment.SequenceAlignment;
import bioinfo.alignment.gotoh.FreeshiftSequenceGotoh;
import bioinfo.alignment.kerbsch.temp.GlobalAngleAligner;
import bioinfo.alignment.kerbsch.temp.InitClass;
import bioinfo.alignment.kerbsch.temp.PairReader;
import bioinfo.alignment.kerbsch.temp.SeqLibrary;


public class aligntest {

	public static void main(String[] args) {
		HashMap<String,char[]> seqlib = SeqLibrary.read("../GoBi_old/referenz/domains.seqlib");
		ArrayList<String[]> pairs = PairReader.parse("vorescore.pairs");
		InitClass init = new InitClass();
		double[][] substMatrix = init.calcGotohInputMatrix(init.calcPolarityScores());
		FreeshiftSequenceGotoh gotoh = new FreeshiftSequenceGotoh(-12.0,-1.0,substMatrix);
		GlobalAngleAligner anglegotoh = new GlobalAngleAligner(-10,-10,"angles","../GoBi_old/DSSP");
		
		Sequence seq1;
		Sequence seq2;
		SequenceAlignment alignment;
//		seq1 = new Sequence("id1","LKAWAAUHREQWMAPQ");
//		seq2 = new Sequence("id2","STACALKJIYCNERY");
//		alignment = anglegotoh.align(seq1, seq2);
		for(String[] pair : pairs){
			seq1 = new Sequence(pair[0],seqlib.get(pair[0]));
			seq2 = new Sequence(pair[1],seqlib.get(pair[1]));
			alignment = anglegotoh.align(seq1, seq2);
			System.out.println(alignment.toStringVerbose());
		}
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
