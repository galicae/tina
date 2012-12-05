package test;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.InputStreamReader;
import java.io.OutputStreamWriter;
import java.util.HashMap;
import java.util.Locale;

import bioinfo.Sequence;
import bioinfo.alignment.SequenceAlignment;
import bioinfo.alignment.gotoh.FreeshiftSequenceGotoh;
import bioinfo.alignment.gotoh.GlobalSequenceGotoh;
import bioinfo.alignment.gotoh.LocalSequenceGotoh;
import bioinfo.alignment.matrices.QuasarMatrix;

public class GotohTest {

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		String seqlibfile = args[1];
		HashMap<String,String> seqlib = new HashMap<String,String>();
		BufferedReader br = null;
		String line = null;
		String[] content;
		String infile = args[0];
		
		try{
			br = new BufferedReader(new InputStreamReader(new FileInputStream(seqlibfile)));
		}catch(Exception e){
			e.printStackTrace();
		}
		if(br == null){
			return;
		}
		try{
			while((line = br.readLine()) != null){
				content = line.split(":");
				seqlib.put(content[0].trim(),content[1].trim());
			}
			br.close();
		}catch(Exception e){
			e.printStackTrace();
		}
		
		

		Sequence seq1;
		Sequence seq2;
		String ali1;
		String ali2;
		SequenceAlignment ali;
		double[][] matrix = QuasarMatrix.parseMatrix(args[2]);
		GlobalSequenceGotoh gotoh = new GlobalSequenceGotoh(-12.0d,-1.0d,matrix);
		try{
			br = new BufferedReader(new InputStreamReader(new FileInputStream(infile)));
		}catch(Exception e){
			e.printStackTrace();
		}
		if(br == null){
			return;
		}
		try{
			while((line = br.readLine()) != null){
				ali1 = br.readLine();
				ali2 = br.readLine();
				line = line.substring(1);
				content = line.split("\\s");
				seq1 = new Sequence(content[0].trim(),seqlib.get(content[0].trim()));
				seq2 = new Sequence(content[1].trim(),seqlib.get(content[1].trim()));
				ali = gotoh.align(seq1, seq2);
				//gotoh.check(ali);
				System.out.println(">"+line+" "+String.format(Locale.US,"%.3f",ali.getScore()));
				System.out.println(ali1);
				System.out.println(ali2);
				System.out.println(content[0].trim()+": "+ali.getRowAsString(0));
				System.out.println(content[1].trim()+": "+ali.getRowAsString(1));
				//gotoh.streamMatricesAsHtml(new BufferedWriter(new OutputStreamWriter(new FileOutputStream("/Users/andreseitz/Desktop/test.ali"))));
			}
			br.close();
		}catch(Exception e){
			e.printStackTrace();
		}
		
	}

}
