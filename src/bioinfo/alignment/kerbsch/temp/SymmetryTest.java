package bioinfo.alignment.kerbsch.temp;

import java.io.BufferedWriter;
import java.io.IOException;
import java.io.OutputStreamWriter;

import bioinfo.Sequence;
import bioinfo.alignment.SequenceAlignment;
import bioinfo.alignment.gotoh.LocalSequenceGotoh;
import bioinfo.alignment.kerbsch.FBGotoh;
import bioinfo.alignment.kerbsch.Kerbsch;
import bioinfo.alignment.matrices.QuasarMatrix;

public class SymmetryTest {

	public static void main(String[] args) throws IOException {
        Sequence seq1 = new Sequence("id1","EEKRNRAITA");
        Sequence seq2 = new Sequence("id2","GERRRSQLDRDQCAYCKEK");
        LocalSequenceGotoh gotoh = new LocalSequenceGotoh(-12.0,-1.0,QuasarMatrix.DAYHOFF_MATRIX);
        BufferedWriter writer = new BufferedWriter(new OutputStreamWriter(System.out));
        InitClass init = new InitClass();
        
        
		Kerbsch kerbsch = new Kerbsch(-1,-1,QuasarMatrix.DAYHOFF_MATRIX,init.calcGotohInputMatrix(init.calcHydropathyScores()),init.calcGotohInputMatrix(init.calcPolarityScores()),SecStructScores.matrix,sec_seqlib);
        FBGotoh fbgotoh = new FBGotoh(-12.0,-1.0,0.38,15.0,63735,QuasarMatrix.DAYHOFF_MATRIX);
        SequenceAlignment normal = gotoh.align(seq1, seq2);
        System.out.println(normal.toStringVerbose());
        fbgotoh.align(seq1.getSequence(), seq2.getSequence());
        writer.close();
	}
        
//        SequenceAlignment normal = gotoh.align(seq1, seq2);
//        System.out.println("normal: "+normal.getScore());
//        System.out.println(normal.getRowAsString(0));
//        System.out.println(normal.getRowAsString(1));
////        gotoh.streamMatricesAsTxt(writer);
//        int[][] normalM = gotoh.getM();
//        
//        SequenceAlignment reverse = gotoh.align(seq1_rev, seq2_rev);
//        System.out.println("reverse: "+reverse.getScore());
//        System.out.println(reverse.getRowAsString(0));
//        System.out.println(reverse.getRowAsString(1));
////        gotoh.streamMatricesAsTxt(writer);
//        int[][] reverseM = gotoh.getM();
//        
//        int[][] hybridM = new int[normalM[0].length][normalM.length];
//        for (int i = 0; i < normalM[0].length; i++) {
//			for (int j = 0; j < normalM.length; j++) {
//				hybridM[i][j] = normalM[i][j] + reverseM[reverseM[0].length-1-i][reverseM.length-1-j];
//			}
//		}
//        
//        for (int i = 0; i < hybridM.length; i++) {
//        	for (int j = 0; j < hybridM.length; j++) {
//				System.out.print("\t"+hybridM[i][j]);
//			}
//        	System.out.println();
//		}
//	}

}
