package bioinfo.alignment.kerbsch;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;


import highscorealignments.CathScopEntry;
import highscorealignments.CathScopHash;
import bioinfo.Sequence;
import bioinfo.alignment.Aligner;
import bioinfo.alignment.SequenceAlignment;
import bioinfo.alignment.gotoh.FreeshiftSequenceGotoh;
import bioinfo.alignment.matrices.QuasarMatrix;

public class AlignmentBenchmarker {

	// different shit
	double go;
	double ge;
	String seqlibfile;
	String pairfile;
	String cathscoplib;
	InitClass matrices = new InitClass();
	double[][] substMatrix;
	HashMap<String, char[]> seqlib;
	ArrayList<String[]> pairs;
	HashMap<String, CathScopEntry> cathscopinfo;
	HashMap<String, Integer> idToIndex = new HashMap<String, Integer>();
	double[][] alignments;
	Aligner gotoh;
	BufferedWriter out = null;

	// family recognition test
	private int cath_recsamefam;
	private int cath_recsamesup;
	private int cath_recsamefold;
	private int cath_recdiffold;
	private int scop_recsamefam;
	private int scop_recsamesup;
	private int scop_recsamefold;
	private int scop_recdiffold;

	// misclassification tests
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
	
	//tempscores for statistic method
	private double maxscore = Double.NEGATIVE_INFINITY;
	private double samefammaxscore_cath = Double.NEGATIVE_INFINITY;
	private double samesupmaxscore_cath = Double.NEGATIVE_INFINITY;
	private double samefoldmaxscore_cath = Double.NEGATIVE_INFINITY;
	private double diffoldmaxscore_cath = Double.NEGATIVE_INFINITY;
	private double samefammaxscore_scop = Double.NEGATIVE_INFINITY;
	private double samesupmaxscore_scop = Double.NEGATIVE_INFINITY;
	private double samefoldmaxscore_scop = Double.NEGATIVE_INFINITY;
	private double diffoldmaxscore_scop = Double.NEGATIVE_INFINITY;

	public AlignmentBenchmarker(String args[], HashMap<String,char[]> sl) {
		go = Double.parseDouble(args[0]);
		ge = Double.parseDouble(args[1]);
		if (args[2].equals("polarity")) {
			substMatrix = matrices.calcGotohInputMatrix(matrices
					.calcPolarityScores());
		} else if (args[2].equals("hydrophob")) {
			substMatrix = matrices.calcGotohInputMatrix(matrices
					.calcHydropathyScores());
		} else if (args[2].equals("secstruct")) {
			substMatrix = SecStructScores.matrix;
		} else if (args[2].equals("sequence")) {
			substMatrix = QuasarMatrix.DAYHOFF_MATRIX;
		}
		gotoh = new FreeshiftSequenceGotoh(go, ge, substMatrix);

		seqlib = sl;
		pairs = PairReader.parse(args[4]);
		cathscopinfo = CathScopHash.read(args[5]);
		alignments = new double[seqlib.size()][seqlib.size()];
		init();
		
		try {
			out = new BufferedWriter(new FileWriter(args[6]));
		} catch (IOException e) {
			System.out.println("cannot initialize writer! (Benchmarker)");
		}	
	}
	
	private void init(){
		for (int i = 0; i < alignments[0].length; i++) {
			for (int j = 0; j < alignments.length; j++) {
				alignments[i][j] = Double.NEGATIVE_INFINITY;
			}
		}
	}
	
	private void statistic(CathScopEntry besthit, CathScopEntry query){
		if(besthit != null){
			System.out.println("found max");
			// for cath family recognition test
			if (query.getCathClazz() == besthit.getCathClazz()
					&& query.getCathFold() == besthit.getCathFold()
					&& query.getCathSupFam() == besthit.getCathSupFam()
					&& query.getCathFam() == besthit.getCathFam()) {
				cath_recsamefam++;
			} else if (query.getCathClazz() == besthit.getCathClazz()
					&& query.getCathFold() == besthit.getCathFold()
					&& query.getCathSupFam() == besthit.getCathSupFam()) {
				cath_recsamesup++;
			} else if (query.getCathClazz() == besthit.getCathClazz()
					&& query.getCathFold() == besthit.getCathFold()) {
				cath_recsamefold++;
			} else {
				cath_recdiffold++;
			}
	
			// for scop family recognition test
			if (query.getScopClazz() == besthit.getScopClazz()
					&& query.getScopFold() == besthit.getScopFold()
					&& query.getScopSupFam() == besthit.getScopSupFam()
					&& query.getScopFam() == besthit.getScopFam()) {
				scop_recsamefam++;
			} else if (query.getScopClazz() == besthit.getScopClazz()
					&& query.getScopFold() == besthit.getScopFold()
					&& query.getScopSupFam() == besthit.getScopSupFam()) {
				scop_recsamesup++;
			} else if (query.getScopClazz() == besthit.getScopClazz()
					&& query.getScopFold() == besthit.getScopFold()) {
				scop_recsamefold++;
			} else {
				scop_recdiffold++;
			}
		}
		
		// cath misclassification test
		if (samefammaxscore_cath >= diffoldmaxscore_cath
				&& diffoldmaxscore_cath != Double.NEGATIVE_INFINITY) {
			cath_missamefam++;
		} else if (samefammaxscore_cath < diffoldmaxscore_cath
				&& samefammaxscore_cath != Double.NEGATIVE_INFINITY) {
			cath_misdiffoldfam++;
		}
		if (samesupmaxscore_cath >= diffoldmaxscore_cath
				&& diffoldmaxscore_cath != Double.NEGATIVE_INFINITY) {
			cath_missamesup++;
		} else if (samesupmaxscore_cath < diffoldmaxscore_cath
				&& samesupmaxscore_cath != Double.NEGATIVE_INFINITY) {
			cath_misdiffoldsup++;
		}
		if (samefoldmaxscore_cath >= diffoldmaxscore_cath
				&& diffoldmaxscore_cath != Double.NEGATIVE_INFINITY) {
			cath_missamefold++;
		} else if (samefoldmaxscore_cath < diffoldmaxscore_cath
				&& samefoldmaxscore_cath != Double.NEGATIVE_INFINITY) {
			cath_misdiffoldfold++;
		}
	
		// scop misclassification test
		if (samefammaxscore_scop >= diffoldmaxscore_scop
				&& diffoldmaxscore_scop != Double.NEGATIVE_INFINITY) {
			scop_missamefam++;
		} else if (samefammaxscore_scop < diffoldmaxscore_scop
				&& samefammaxscore_scop != Double.NEGATIVE_INFINITY) {
			scop_misdiffoldfam++;
		}
		if (samesupmaxscore_scop >= diffoldmaxscore_scop
				&& diffoldmaxscore_scop != Double.NEGATIVE_INFINITY) {
			scop_missamesup++;
		} else if (samesupmaxscore_scop < diffoldmaxscore_scop
				&& samesupmaxscore_scop != Double.NEGATIVE_INFINITY) {
			scop_misdiffoldsup++;
		}
		if (samefoldmaxscore_scop >= diffoldmaxscore_scop
				&& diffoldmaxscore_scop != Double.NEGATIVE_INFINITY) {
			scop_missamefold++;
		} else if (samefoldmaxscore_scop < diffoldmaxscore_scop
				&& samefoldmaxscore_scop != Double.NEGATIVE_INFINITY) {
			scop_misdiffoldfold++;
		}
	}
	
	private void resetStatistic(){
		maxscore = Double.NEGATIVE_INFINITY;
		samefammaxscore_cath = Double.NEGATIVE_INFINITY;
		samesupmaxscore_cath = Double.NEGATIVE_INFINITY;
		samefoldmaxscore_cath = Double.NEGATIVE_INFINITY;
		diffoldmaxscore_cath = Double.NEGATIVE_INFINITY;
		samefammaxscore_scop = Double.NEGATIVE_INFINITY;
		samesupmaxscore_scop = Double.NEGATIVE_INFINITY;
		samefoldmaxscore_scop = Double.NEGATIVE_INFINITY;
		diffoldmaxscore_scop = Double.NEGATIVE_INFINITY;
	}

	public void benchmark() {
		double score = Double.NEGATIVE_INFINITY;
		resetStatistic();
		
		String lastquery = "";
		SequenceAlignment temp = null;
		CathScopEntry query = null;
		CathScopEntry template;
		CathScopEntry besthit = null;

		int index = 0;
		int id1;
		int id2;
		
		for (String[] pair : pairs) {

			// give index to ID for storing the score of the alignment in a
			// matrix
			if (!idToIndex.containsKey(pair[0])) {
				idToIndex.put(pair[0], index);
				index++;
			}
			if (!idToIndex.containsKey(pair[1])) {
				idToIndex.put(pair[1], index);
				index++;
			}

			// here comes a new target to thread
			if (!pair[0].equals(lastquery)) {
				statistic(besthit,query);
				resetStatistic();
				besthit = null;
				lastquery = pair[0];
				query = cathscopinfo.get(pair[0]);
			}
			template = cathscopinfo.get(pair[1]);
			
			// make Alignment and check for besthit
			id1 = idToIndex.get(pair[0]);
			id2 = idToIndex.get(pair[1]);
			if (alignments[id1][id2] == Double.NEGATIVE_INFINITY) {
				temp = (SequenceAlignment) gotoh.align(new Sequence(pair[0],
						seqlib.get(pair[0])),
						new Sequence(pair[1], seqlib.get(pair[1])));
				score = temp.getScore();
				alignments[id1][id2] = score;
				alignments[id2][id1] = score;
			}
			//else get score from alignment matrix
			else{
				score = alignments[id1][id2];
			}
			if (score > maxscore) {
				maxscore = score;
				besthit = cathscopinfo.get(pair[1]);
			}
			
			// for cath misclassification tests
			if (query.getCathClazz() == template.getCathClazz()
					&& query.getCathFold() == template.getCathFold()
					&& query.getCathSupFam() == template.getCathSupFam()
					&& query.getCathFam() == template.getCathFam()) {
				if (score > samefammaxscore_cath) {
					samefammaxscore_cath = score;
				}
			}
			if (query.getCathClazz() == template.getCathClazz()
					&& query.getCathFold() == template.getCathFold()
					&& query.getCathSupFam() == template.getCathSupFam()
					&& query.getCathFam() != template.getCathFam()) {
				if (score > samesupmaxscore_cath) {
					samesupmaxscore_cath = score;
				}
			}
			if (query.getCathClazz() == template.getCathClazz()
					&& query.getCathFold() == template.getCathFold()
					&& query.getCathSupFam() != template.getCathSupFam()) {
				if (score > samefoldmaxscore_cath) {
					samefoldmaxscore_cath = score;
				}
			}
			if (query.getCathClazz() != template.getCathClazz()
					|| query.getCathFold() != template.getCathFold()) {
				if (score > diffoldmaxscore_cath) {
					diffoldmaxscore_cath = score;
				}
			}

			// for scop misclassification tests
			if (query.getScopClazz() == template.getScopClazz()
					&& query.getScopFold() == template.getScopFold()
					&& query.getScopSupFam() == template.getScopSupFam()
					&& query.getScopFam() == template.getScopFam()) {
				if (score > samefammaxscore_scop) {
					samefammaxscore_scop = score;
				}
			}
			if (query.getScopClazz() == template.getScopClazz()
					&& query.getScopFold() == template.getScopFold()
					&& query.getScopSupFam() == template.getScopSupFam()
					&& query.getScopFam() != template.getScopFam()) {
				if (score > samesupmaxscore_scop) {
					samesupmaxscore_scop = score;
				}
			}
			if (query.getScopClazz() == template.getScopClazz()
					&& query.getScopFold() == template.getScopFold()
					&& query.getScopSupFam() != template.getScopSupFam()) {
				if (score > samefoldmaxscore_scop) {
					samefoldmaxscore_scop = score;
				}
			}
			if (query.getScopClazz() != template.getScopClazz()
					|| query.getScopFold() != template.getScopFold()) {
				if (score > diffoldmaxscore_scop) {
					diffoldmaxscore_scop = score;
				}
			}
		}
		statistic(besthit,query);
	}
	
//	private void printResultList(CathScopEntry besthit, CathScopEntry query) throws IOException {
//		out.append(query.getID()+"\t"+query.getCathClazz()+"."+query.getCathFold()+"."+query.getCathSupFam()+"."+query.getCathFam()+"\t");
//		out.append(query.getScopClazz()+"."+query.getScopFold()+"."+query.getScopSupFam()+"."+query.getScopFam()+"\t");
//		out.append(besthit.getID()+"\t"+besthit.getCathClazz()+"."+besthit.getCathFold()+"."+besthit.getCathSupFam()+"."+besthit.getCathFam()+"\t");
//		out.append(besthit.getScopClazz()+"."+besthit.getScopFold()+"."+besthit.getScopSupFam()+"."+besthit.getScopFam()+"\t"+alignments[idToIndex.get(query.getID())][idToIndex.get(besthit.getID())]+"\n");
//	}

	public void printResults() throws IOException {
		out.write(cath_recsamefam + "\n" + cath_recsamesup + "\n"
				+ cath_recsamefold + "\n" + cath_recdiffold + "\n");
		// family misclassification
		out.write(cath_missamefam + "\n" + cath_misdiffoldfam + "\n");
		// sup misclassification
		out.write(cath_missamesup + "\n" + cath_misdiffoldsup + "\n");
		// fold misclassification
		out.write(cath_missamefold + "\n" + cath_misdiffoldfold + "\n");

		// scop
		// family recognition
		out.write(scop_recsamefam + "\n" + scop_recsamesup + "\n"
				+ scop_recsamefold + "\n" + scop_recdiffold + "\n");
		// family misclassification
		out.write(scop_missamefam + "\n" + scop_misdiffoldfam + "\n");
		// sup misclassification
		out.write(scop_missamesup + "\n" + scop_misdiffoldsup + "\n");
		// fold misclassification
		out.write(scop_missamefold + "\n" + scop_misdiffoldfold + "\n");
		out.close();
	}
}
