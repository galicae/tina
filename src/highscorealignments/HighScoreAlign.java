package highscorealignments;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map.Entry;

import bioinfo.Sequence;
import bioinfo.alignment.Aligner;
import bioinfo.alignment.SequenceAlignment;

import bioinfo.alignment.gotoh.FreeshiftSequenceGotoh;
import bioinfo.alignment.gotoh.GlobalSequenceGotoh;
import bioinfo.alignment.gotoh.LocalSequenceGotoh;
import bioinfo.alignment.kerbsch.temp.SeqLibrary;
import bioinfo.alignment.matrices.QuasarMatrix;

public class HighScoreAlign {
	
	private HashMap<String,SequenceAlignment> hsalign = new HashMap<String,SequenceAlignment>();
	private HashMap<String,ArrayList<String>> pairs;
	private HashMap<String, char[]> seqlib;
	private double[][] scoringmatrix;
	private double go;
	private double ge;
	private BufferedWriter out;
	
	public HighScoreAlign(String inpairs, String outpairs, String seqlibpath){
		this.pairs = HashPairReader.readPairs(inpairs, outpairs);
		this.seqlib = SeqLibrary.read(seqlibpath);
	}
	
	public void calculation(String matrixpath, int go, int ge, String mode, String outpath){
		this.scoringmatrix = QuasarMatrix.parseMatrix(matrixpath);
		this.go = (double)go;
		this.ge = (double)ge;
		Aligner gotoh;

		try {
			this.out = new BufferedWriter(new FileWriter(outpath));
		} catch (IOException e) {
			System.out.println("cannot make new alignwriter");
		}
		
		if(mode.equals("global")){
			gotoh = new GlobalSequenceGotoh(this.go,this.ge,this.scoringmatrix);
		}
		else if(mode.equals("local")){
			gotoh = new LocalSequenceGotoh(this.go,this.ge,this.scoringmatrix);
		}
		else{
			gotoh = new FreeshiftSequenceGotoh(this.go,this.ge,this.scoringmatrix);
		}
				
		double maxscore;
		SequenceAlignment maxalign;
		SequenceAlignment temp;
		for(Entry<String,ArrayList<String>> actid : this.pairs.entrySet()){
			maxscore = Double.MIN_VALUE;
			maxalign = null;
			for(String alignpartner : actid.getValue()){
				temp = (SequenceAlignment) gotoh.align(new Sequence(actid.getKey(),seqlib.get(actid.getKey())),new Sequence(alignpartner,seqlib.get(alignpartner)));
				if (temp.getScore() > maxscore){
					maxalign = temp;
					maxscore = temp.getScore();
				}
			}
			writeMaxAlign(maxalign);
		}
		try {
			out.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
	
	public HashMap<String,SequenceAlignment> getHighScoreAligns(){
		return this.hsalign;
	}
	
	public void writeMaxAlign(SequenceAlignment align){
		try {
			String id1 = align.getComponent(0).getId();
			String id2 = align.getComponent(1).getId();
			double score = align.getScore();
			out.write(">"+id1+" "+id2+" "+score+"\n");
			out.write(id1+": "+align.getRowAsString(0)+"\n");
			out.write(id2+": "+align.getRowAsString(1)+"\n");
		} catch (IOException e) {
			System.out.println("Cannot write alignments");
		}
	}
//	public static void main (String[] args){
//		HighScoreAlign test = new HighScoreAlign(args[0],args[1],args[2]);
//		test.calculation(args[3], -12, -1, "freeshift","nikohighscore.align");
//	}
}
