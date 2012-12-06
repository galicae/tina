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
import bioinfo.alignment.matrices.QuasarMatrix;

public class BMCathAndScop {
	
	private HashMap<String,char[]> seqlib;
	private HashMap<String,CathScopEntry> cathscopinfos;
	private HashMap<String,ArrayList<String>> pairs;
	
	//family recognition test
	private int cath_recsamefam;
	private int cath_recsamesup;
	private int cath_recsamefold;
	private int cath_recdiffold;
	private int scop_recsamefam;
	private int scop_recsamesup;
	private int scop_recsamefold;
	private int scop_recdiffold;
	
	//misclassification tests
	private int cath_missamefam;
	private int cath_missamesup;
	private int cath_missamefold;
	private int cath_misdiffoldfam;
	private int cath_misdiffoldsup;
	private int cath_misdiffoldfold;
	private int scop_missamefam;
	private int scop_missamesup;
	private int scop_missamefold;
	private int scop_misdiffoldfam;
	private int scop_misdiffoldsup;
	private int scop_misdiffoldfold;
	
	private Aligner gotoh;
	
	private BufferedWriter out = null; 
	
	public BMCathAndScop(String seqlibpath, String cathscopinfopath, String pairlistpath, String matrixpath,int go,int ge, String mode, String outpath){
		try {
			this.out = new BufferedWriter(new FileWriter(outpath));
		} catch (IOException e) {
			System.out.println("geht nischt");
		}
		this.seqlib = SeqLibrary.parse(seqlibpath);
		this.cathscopinfos = CathScopHash.read(cathscopinfopath);
		this.pairs = HashPairReader.readPairs(pairlistpath+"/cathscop.inpairs", pairlistpath+"/cathscop.outpairs");
		double[][] scoringmatrix = QuasarMatrix.parseMatrix(matrixpath);
		
		this.cath_recdiffold = 0;
		this.cath_recsamefam = 0;
		this.cath_recsamefold = 0;
		this.cath_recsamesup = 0;
		this.scop_recdiffold = 0;
		this.scop_recsamefam = 0;
		this.scop_recsamefold = 0;
		this.scop_recsamesup = 0;
		
		this.cath_misdiffoldfam = 0;
		this.cath_misdiffoldsup = 0;
		this.cath_misdiffoldfold = 0;
		this.cath_missamefam = 0;
		this.cath_missamesup = 0;
		this.cath_missamefold = 0;
		this.scop_misdiffoldfam = 0;
		this.scop_misdiffoldsup = 0;
		this.scop_misdiffoldfold = 0;
		this.scop_missamefam = 0;
		this.scop_missamesup = 0;
		this.scop_missamefold = 0;

		if(mode.equals("global")){
			gotoh = new GlobalSequenceGotoh(go,ge,scoringmatrix);
		}
		else if(mode.equals("local")){
			gotoh = new LocalSequenceGotoh(go,ge,scoringmatrix);
		}
		else{
			gotoh = new FreeshiftSequenceGotoh(go,ge,scoringmatrix);
		}
		benchmark();
	}
	
	public void benchmark(){
		double maxscore;
		double score;
		double samefammaxscore_cath;
		double samesupmaxscore_cath;
		double samefoldmaxscore_cath;
		double diffoldmaxscore_cath;
		double samefammaxscore_scop;
		double samesupmaxscore_scop;
		double samefoldmaxscore_scop;
		double diffoldmaxscore_scop;
		
		for(Entry<String,ArrayList<String>> query_entry : pairs.entrySet()){
		CathScopEntry besthit = null;
		maxscore = Double.MIN_VALUE;
		samefammaxscore_cath = Double.MIN_VALUE;
		samesupmaxscore_cath = Double.MIN_VALUE;
		samefoldmaxscore_cath = Double.MIN_VALUE;
		diffoldmaxscore_cath = Double.MIN_VALUE;
		samefammaxscore_scop = Double.MIN_VALUE;
		samesupmaxscore_scop = Double.MIN_VALUE;
		samefoldmaxscore_scop = Double.MIN_VALUE;
		diffoldmaxscore_scop = Double.MIN_VALUE;
		SequenceAlignment temp;
		CathScopEntry query = cathscopinfos.get(query_entry.getKey());
		
			for(String template : query_entry.getValue()){
				temp = (SequenceAlignment)gotoh.align(new Sequence(query.getID(),seqlib.get(query.getID())),new Sequence(template,seqlib.get(template)));
				score = temp.getScore();
				if(score > maxscore){
					maxscore = score;
					besthit = cathscopinfos.get(template);
				}
				//for cath misclassification tests
				if(query.getCathClazz() == besthit.getCathClazz() && query.getCathFold() == besthit.getCathFold() && query.getCathSupFam() == besthit.getCathSupFam() && query.getCathFam() == besthit.getCathFam()){
					if(score > samefammaxscore_cath){
						samefammaxscore_cath = score;
					}
				}
				if(query.getCathClazz() == besthit.getCathClazz() && query.getCathFold() == besthit.getCathFold() && query.getCathSupFam() == besthit.getCathSupFam()){
					if(score > samesupmaxscore_cath){
						samesupmaxscore_cath = score;
					}
				}
				if(query.getCathClazz() == besthit.getCathClazz() && query.getCathFold() == besthit.getCathFold()){
					if(score > samefoldmaxscore_cath){
						samefoldmaxscore_cath = score;
					}
				}
				if(query.getCathClazz() != besthit.getCathClazz() || query.getCathFold() != besthit.getCathFold()){
					if(score > diffoldmaxscore_cath){
						diffoldmaxscore_cath = score;
					}
				}

				//for scop misclassification tests
				if(query.getScopClazz() == besthit.getScopClazz() && query.getScopFold() == besthit.getScopFold() && query.getScopSupFam() == besthit.getScopSupFam() && query.getScopFam() == besthit.getScopFam()){
					if(score > samefammaxscore_scop){
						samefammaxscore_scop = score;
					}
				}
				if(query.getScopClazz() == besthit.getScopClazz() && query.getScopFold() == besthit.getScopFold() && query.getScopSupFam() == besthit.getScopSupFam()){
					if(score > samesupmaxscore_scop){
						samesupmaxscore_scop = score;
					}
				}
				if(query.getScopClazz() == besthit.getScopClazz() && query.getScopFold() == besthit.getScopFold()){
					if(score > samefoldmaxscore_scop){
						samefoldmaxscore_scop = score;
					}
				}
				if(query.getScopClazz() != besthit.getScopClazz() || query.getScopFold() != besthit.getScopFold()){
					if(score > diffoldmaxscore_scop){
						diffoldmaxscore_scop = score;
					}
				}
			}
			
			//for cath family recognition test
			if(query.getCathClazz() == besthit.getCathClazz() && query.getCathFold() == besthit.getCathFold() && query.getCathSupFam() == besthit.getCathSupFam() && query.getCathFam() == besthit.getCathFam()){
				cath_recsamefam++;
			}
			else if(query.getCathClazz() == besthit.getCathClazz() && query.getCathFold() == besthit.getCathFold() && query.getCathSupFam() == besthit.getCathSupFam()){
				cath_recsamesup++;
			}
			else if(query.getCathClazz() == besthit.getCathClazz() && query.getCathFold() == besthit.getCathFold()){
				cath_recsamefold++;
			}
			else{
				cath_recdiffold++;
			}
			//for scop family recognition test
			if(query.getScopClazz() == besthit.getScopClazz() && query.getScopFold() == besthit.getScopFold() && query.getScopSupFam() == besthit.getScopSupFam() && query.getScopFam() == besthit.getScopFam()){
				scop_recsamefam++;
			}
			else if(query.getScopClazz() == besthit.getScopClazz() && query.getScopFold() == besthit.getScopFold() && query.getScopSupFam() == besthit.getScopSupFam()){
				scop_recsamesup++;
			}
			else if(query.getScopClazz() == besthit.getScopClazz() && query.getScopFold() == besthit.getScopFold()){
				scop_recsamefold++;
			}
			else{
				scop_recdiffold++;
			}
			
			//cath misclassification test
			if(samefammaxscore_cath > diffoldmaxscore_cath){
				cath_missamefam++;
			} else {cath_misdiffoldfam++;}
			if(samesupmaxscore_cath > diffoldmaxscore_cath){
				cath_missamesup++;
			} else {cath_misdiffoldsup++;}
			if(samefoldmaxscore_cath > diffoldmaxscore_cath){
				cath_missamefold++;
			} else {cath_misdiffoldfold++;}
			//scop misclassification test
			if(samefammaxscore_scop > diffoldmaxscore_scop){
				scop_missamefam++;
			} else {scop_misdiffoldfam++;}
			if(samesupmaxscore_scop > diffoldmaxscore_scop){
				scop_missamesup++;
			} else {scop_misdiffoldsup++;}
			if(samefoldmaxscore_scop > diffoldmaxscore_scop){
				scop_missamefold++;
			} else {scop_misdiffoldfold++;}
		}
		
		try {
			printResults();
			out.close();
		} catch (IOException e) {
			System.out.println("Cannot write results");
		}
	}
	
	public void printResults() throws IOException{
		out.write(cath_recsamefam+"\n"+cath_recsamesup+"\n"+cath_recsamefold+"\n"+cath_recdiffold);
		//family misclassification
		out.write(cath_missamefam+"\n"+cath_misdiffoldfam);
		//sup misclassification
		out.write(cath_missamesup+"\n"+cath_misdiffoldsup);
		//fold misclassification
		out.write(cath_missamefold+"\n"+cath_misdiffoldfold);

		//scop
		//family recognition
		out.write(scop_recsamefam+"\n"+scop_recsamesup+"\n"+scop_recsamefold+"\n"+scop_recdiffold);
		//family misclassification
		out.write(scop_missamefam+"\n"+scop_misdiffoldfam);
		//sup misclassification
		out.write(scop_missamesup+"\n"+scop_misdiffoldsup);
		//fold misclassification
		out.write(scop_missamefold+"\n"+scop_misdiffoldfold);	
	}
	
	public static void main(String[] args){
		BMCathAndScop benchmark = new BMCathAndScop(args[0], args[1], args[2], args[3], Integer.parseInt(args[4]), Integer.parseInt(args[5]), args[6], args[7]);
	}
}
