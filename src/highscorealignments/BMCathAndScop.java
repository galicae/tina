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

public class BMCathAndScop {
	
	private HashMap<String,char[]> seqlib;
	private HashMap<String,CathScopEntry> cathinfos;
	private HashMap<String,ArrayList<String>> pairs;
	
	//family recognition test
	private int recsamefam;
	private int recsamesup;
	private int recsamefold;
	private int recdiffold;
	
	//misclassification tests
	private int missamefam;
	private int missamesup;
	private int missamefold;
	private int misdiffoldfam;
	private int misdiffoldsup;
	private int misdiffoldfold;
	
	private Aligner gotoh;
	
	public BMCathAndScop(HashMap<String,char[]> seqlib, HashMap<String,CathScopEntry> cathinfos, HashMap<String,ArrayList<String>> pairlist, String matrixpath,int go,int ge, String mode){
		this.seqlib = seqlib;
		this.cathinfos = cathinfos;
		this.recdiffold = 0;
		this.recsamefam = 0;
		this.recsamefold = 0;
		this.recsamesup = 0;
		
		this.misdiffoldfam = 0;
		this.misdiffoldsup = 0;
		this.misdiffoldfold = 0;
		this.missamefam = 0;
		this.missamesup = 0;
		this.missamefold = 0;

		this.pairs = pairlist;
		double[][] scoringmatrix = QuasarMatrix.parseMatrix(matrixpath);
		
		if(mode.equals("global")){
			gotoh = new GlobalSequenceGotoh(go,ge,scoringmatrix);
		}
		else if(mode.equals("local")){
			gotoh = new LocalSequenceGotoh(go,ge,scoringmatrix);
		}
		else{
			gotoh = new FreeshiftSequenceGotoh(go,ge,scoringmatrix);
		}
	}
	
	public void benchmark(){
		double maxscore;
		double score;
		double samefammaxscore;
		double samesupmaxscore;
		double samefoldmaxscore;
		double diffoldmaxscore;
		
		for(Entry<String,ArrayList<String>> query_entry : pairs.entrySet()){
		CathScopEntry besthit = null;
		maxscore = Double.MIN_VALUE;
		samefammaxscore = Double.MIN_VALUE;
		samesupmaxscore = Double.MIN_VALUE;
		samefoldmaxscore = Double.MIN_VALUE;
		diffoldmaxscore = Double.MIN_VALUE;
		SequenceAlignment temp;
		CathScopEntry query = cathinfos.get(query_entry.getKey());
		
			for(String template : query_entry.getValue()){
				temp = (SequenceAlignment)gotoh.align(new Sequence(query.getID(),seqlib.get(query.getID())),new Sequence(template,seqlib.get(template)));
				score = temp.getScore();
				if(score > maxscore){
					maxscore = score;
					besthit = cathinfos.get(template);
				}
				
				//for misclassification tests
				if(query.getCathClazz() == besthit.getCathClazz() && query.getCathFold() == besthit.getCathFold() && query.getCathSupFam() == besthit.getCathSupFam() && query.getCathFam() == besthit.getCathFam()){
					if(score > samefammaxscore){
						samefammaxscore = score;
					}
				}
				else if(query.getCathClazz() == besthit.getCathClazz() && query.getCathFold() == besthit.getCathFold() && query.getCathSupFam() == besthit.getCathSupFam()){
					if(score > samesupmaxscore){
						samesupmaxscore = score;
					}
				}
				else if(query.getCathClazz() == besthit.getCathClazz() && query.getCathFold() == besthit.getCathFold()){
					if(score > samefoldmaxscore){
						samefoldmaxscore = score;
					}
				}
				else{
					if(score > diffoldmaxscore){
						diffoldmaxscore = score;
					}
				}			
			}
			
			//for family recognition test
			if(query.getCathClazz() == besthit.getCathClazz() && query.getCathFold() == besthit.getCathFold() && query.getCathSupFam() == besthit.getCathSupFam() && query.getCathFam() == besthit.getCathFam()){
				recsamefam++;
			}
			else if(query.getCathClazz() == besthit.getCathClazz() && query.getCathFold() == besthit.getCathFold() && query.getCathSupFam() == besthit.getCathSupFam()){
				recsamesup++;
			}
			else if(query.getCathClazz() == besthit.getCathClazz() && query.getCathFold() == besthit.getCathFold()){
				recsamefold++;
			}
			else{
				recdiffold++;
			}
			
			//misclassification test
			if(samefammaxscore > diffoldmaxscore){
				missamefam++;
			} else {misdiffoldfam++;}
			if(samesupmaxscore > diffoldmaxscore){
				missamesup++;
			} else {misdiffoldsup++;}
			if(samefoldmaxscore > diffoldmaxscore){
				missamefold++;
			} else {misdiffoldfold++;}
		}
	}
	
	public void printResults(){
		//family recognition
		System.out.println(recsamefam+"\n"+recsamesup+"\n"+recsamefold+"\n"+recdiffold);
		//family misclassification
		System.out.println(missamefam+"\n"+misdiffoldfam);
		//sup misclassification
		System.out.println(missamesup+"\n"+misdiffoldsup);
		//fold misclassification
		System.out.println(missamefold+"\n"+misdiffoldfold);

	}
}
