package test;


import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.HashMap;

import bioinfo.Sequence;
import bioinfo.alignment.SequenceAlignment;
import bioinfo.alignment.d123.FreeshiftSequence123D;
import bioinfo.alignment.gotoh.Gotoh;
import bioinfo.alignment.kerbsch.temp.SeqLibrary;
import bioinfo.alignment.matrices.MatrixReader123D;
import bioinfo.alignment.matrices.QuasarMatrix;
import bioinfo.proteins.SSCCEntry;
import bioinfo.proteins.SSCCReader;

public class OneTwoThreeDTest {

	public static void main(String[] args) {
		
		//initialize
		BufferedReader in;
		HashMap<String,char[]> seqlib = SeqLibrary.read(args[7]);
		double[][] scoringmatrix = QuasarMatrix.parseMatrix(args[0]);
		double[][][] potentials = MatrixReader123D.readPotentials(args[1], args[2], args[3]);
		double[][] secstructpref = MatrixReader123D.readSecStructPref(args[4]);
		double[][] weights = MatrixReader123D.readWeights(args[5]);
		String ssccpath = args[6];
		FreeshiftSequence123D threader = new FreeshiftSequence123D(-12.0,-1.0,scoringmatrix,secstructpref,weights,potentials);
		
		//run that shit
		try{
			in = new BufferedReader(new FileReader(args[8]));
			String line;
			String[] temp;
			Sequence seq1;
			Sequence seq2;
			SequenceAlignment seqalign;
			SSCCEntry ssccentry;
			while((line = in.readLine()) != null){
				temp = line.split("\\s+");	
				ssccentry = SSCCReader.readSSCC(ssccpath+"/"+temp[1]+".sscc");
				seq1 = new Sequence(temp[0],seqlib.get(temp[0]));
				seq2 = new Sequence(temp[1],seqlib.get(temp[1]));
				seqalign = threader.align(seq1, seq2,ssccentry);
				threader.streamMatricesAsHtml(new BufferedWriter(new FileWriter("matrices.html")),Gotoh.FACTOR*Gotoh.FACTOR);
				System.out.println(temp[0]+"\t"+temp[1]+"\t"+seqalign.getScore());
				System.out.println(temp[0]+": "+seqalign.getRowAsString(0));
				System.out.println(temp[1]+": "+seqalign.getRowAsString(1));
			}
			in.close();
		} catch(IOException e){
			System.out.println("cannot read pair-file");
		}
		
	}

}
