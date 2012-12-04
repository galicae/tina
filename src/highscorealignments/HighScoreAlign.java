package highscorealignments;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map.Entry;

import bioinfo.Sequence;
import bioinfo.alignment.Aligner;
import bioinfo.alignment.FreeshiftSequenceGotoh;
import bioinfo.alignment.GlobalSequenceGotoh;
import bioinfo.alignment.LocalSequenceGotoh;
import bioinfo.alignment.SequenceAlignment;

import bioinfo.alignment.matrices.QuasarMatrix;

public class HighScoreAlign {
	
	private HashMap<String,SequenceAlignment> hsalign = new HashMap<String,SequenceAlignment>();
	private ArrayList<CathScopEntry[]> inpairs;
	private HashMap<String,ArrayList<CathScopEntry[]>> outpairs;
	private HashMap<String, char[]> seqlib;
	private double[][] scoringmatrix;
	private double go;
	private double ge;
	
	public HighScoreAlign(String inpairs, String outpairs, String seqlibpath){
		this.inpairs = PairReaderCathScop.readPairs(inpairs);
		this.outpairs = HashPairReaderCathScop.readPairs(outpairs);
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
		
		String lastcath = "";
		String actcath;
		double maxscore = Double.MIN_VALUE;
		SequenceAlignment temp;
		for(CathScopEntry[] pair : inpairs){
			actcath = pair[0].getID();
			temp = (SequenceAlignment) gotoh.align(new Sequence(pair[0].getID(),seqlib.get(pair[0].getID())),new Sequence(pair[1].getID(),seqlib.get(pair[1].getID())));
			
			if(lastcath.equals(actcath)){
				if (temp.getScore() > maxscore){
					hsalign.put(actcath,(SequenceAlignment) gotoh.align(new Sequence(pair[0].getID(),seqlib.get(pair[0].getID())),new Sequence(pair[1].getID(),seqlib.get(pair[1].getID()))));
					maxscore = temp.getScore();
				}
			}
			else{
				lastcath = actcath;
				hsalign.put(actcath,(SequenceAlignment) gotoh.align(new Sequence(pair[0].getID(),seqlib.get(pair[0].getID())),new Sequence(pair[1].getID(),seqlib.get(pair[1].getID()))));
				maxscore = hsalign.get(actcath).getScore();
				if(this.outpairs.containsKey(actcath)){
					for(CathScopEntry[] out : this.outpairs.get(actcath)){
						SequenceAlignment temp2 = (SequenceAlignment) gotoh.align(new Sequence(out[0].getID(),this.seqlib.get(out[0].getID())), new Sequence(out[1].getID(),this.seqlib.get(out[1].getID())));
						if(temp2.getScore() > maxscore){
							maxscore = temp2.getScore();
							hsalign.put(actcath,(SequenceAlignment) gotoh.align(new Sequence(out[0].getID(),this.seqlib.get(out[0].getID())), new Sequence(out[1].getID(),this.seqlib.get(out[1].getID()))));
						}
					}
				}
			}		
		}
	}
	
	public HashMap<String,SequenceAlignment> getHighScoreAligns(){
		return this.hsalign;
	}
	
	public static void main (String[] args){
		HighScoreAlign test = new HighScoreAlign(args[0],args[1],args[2]);
		test.calculation(args[3], -12, -1, "freeshift");
		for(Entry<String,SequenceAlignment> a : test.getHighScoreAligns().entrySet()){
			System.out.println(a.getValue().toStringVerbose());
		}
	}
}
