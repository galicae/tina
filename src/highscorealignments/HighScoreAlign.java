package highscorealignments;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map.Entry;

import bioinfo.Sequence;
import bioinfo.alignment.Aligner;
import bioinfo.alignment.SequenceAlignment;

import bioinfo.alignment.gotoh.FreeshiftSequenceGotoh;
import bioinfo.alignment.gotoh.GlobalSequenceGotoh;
import bioinfo.alignment.gotoh.LocalSequenceGotoh;
import bioinfo.alignment.matrices.QuasarMatrix;

public class HighScoreAlign {
	
	private HashMap<String,SequenceAlignment> hsalign = new HashMap<String,SequenceAlignment>();
	private HashMap<String,ArrayList<String>> pairs;
	private HashMap<String, char[]> seqlib;
	private double[][] scoringmatrix;
	private double go;
	private double ge;
	
	public HighScoreAlign(String inpairs, String outpairs, String seqlibpath){
		this.pairs = HashPairReader.readPairs(inpairs, outpairs);
		this.seqlib = SeqLibrary.parse(seqlibpath);	
	}
	
	public void calculation(String matrixpath, int go, int ge, String mode){
		this.scoringmatrix = QuasarMatrix.parseMatrix(matrixpath);
		this.go = (double)go;
		this.ge = (double)ge;
		Aligner gotoh;
		if(mode.equals("global")){
			gotoh = new GlobalSequenceGotoh(this.go,this.ge,this.scoringmatrix);
		}
		else if(mode.equals("local")){
			gotoh = new LocalSequenceGotoh(this.go,this.ge,this.scoringmatrix);
		}
		else{
			gotoh = new FreeshiftSequenceGotoh(this.go,this.ge,this.scoringmatrix);
		}
				
		double maxscore = Double.MIN_VALUE;
		SequenceAlignment temp;
		for(Entry<String,ArrayList<String>> actid : this.pairs.entrySet()){
			for(String alignpartner : actid.getValue()){
				temp = (SequenceAlignment) gotoh.align(new Sequence(actid.getKey(),seqlib.get(actid.getKey())),new Sequence(alignpartner,seqlib.get(alignpartner)));
				if (temp.getScore() > maxscore){
					hsalign.put(actid.getKey(),temp.duplicate());
					maxscore = temp.getScore();
				}
			}
		}
	}
	
	public HashMap<String,SequenceAlignment> getHighScoreAligns(){
		return this.hsalign;
	}
	
//	public static void main (String[] args){
//		HighScoreAlign test = new HighScoreAlign(args[0],args[1],args[2]);
//		test.calculation(args[3], -12, -1, "freeshift");
//		for(Entry<String,SequenceAlignment> a : test.getHighScoreAligns().entrySet()){
//			System.out.println(a.getValue().toStringVerbose());
//		}
//	}
}
