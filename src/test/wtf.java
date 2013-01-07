package test;

import bioinfo.Sequence;
import bioinfo.alignment.SequenceAlignment;
import bioinfo.alignment.gotoh.FreeshiftSequenceGotoh;
import bioinfo.alignment.matrices.QuasarMatrix;

public class wtf {


	public static void main(String[] args) {
		Sequence seq1 = new Sequence("id1","ABCCDEF");
		Sequence seq2 = new Sequence("id2","CCABBFE");
		Sequence seq1_rev = new Sequence("id3","FEDCCBA");
		Sequence seq2_rev = new Sequence("id4","EFBBACC");
		FreeshiftSequenceGotoh gotoh = new FreeshiftSequenceGotoh(-12.0,-1.0,QuasarMatrix.DAYHOFF_MATRIX);
		SequenceAlignment normal = gotoh.align(seq1, seq2);
		SequenceAlignment reverse = gotoh.align(seq1_rev, seq2_rev);
		
		System.out.println("normal: "+normal.getScore());
		System.out.println("reverse: "+reverse.getScore());
	}

}
