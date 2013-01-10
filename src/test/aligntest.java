package test;

import java.util.ArrayList;
import java.util.HashMap;

import bioinfo.Sequence;
import bioinfo.alignment.SequenceAlignment;
import bioinfo.alignment.gotoh.FreeshiftSequenceGotoh;
import bioinfo.alignment.kerbsch.InitClass;

import highscorealignments.SeqLibrary;

public class aligntest {

	public static void main(String[] args) {
		HashMap<String,char[]> seqlib = SeqLibrary.parse("../Gobi_old/referenz/domains.seqlib");
		ArrayList<String[]> pairs = PairReader.parse("benchmark.pairs");
		InitClass init = new InitClass();
		double[][] substMatrix = init.calcGotohInputMatrix(init.calcPolarityScores());
		FreeshiftSequenceGotoh gotoh = new FreeshiftSequenceGotoh(-12.0,-1.0,substMatrix);
		
		Sequence seq1;
		Sequence seq2;
		SequenceAlignment alignment;
		for(String[] pair : pairs){
			seq1 = new Sequence(pair[0],seqlib.get(pair[0]));
			seq2 = new Sequence(pair[1],seqlib.get(pair[1]));
			alignment = gotoh.align(seq1, seq2);
			System.out.println(alignment.toStringVerbose());
		}
	}

}