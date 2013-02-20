package bioinfo.alignment.kerbsch;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.Reader;
import java.util.ArrayList;
import java.util.HashMap;

import org.jmol.minimize.forcefield.ForceField;

import bioinfo.Sequence;
import bioinfo.alignment.Aligner;
import bioinfo.alignment.SequenceAlignment;
import bioinfo.alignment.kerbsch.temp.InitClass;
import bioinfo.alignment.kerbsch.temp.SecStructScores;
import bioinfo.alignment.kerbsch.temp.SeqLibrary;
import bioinfo.alignment.matrices.QuasarMatrix;

public class ValidateAlignments {

	// different shit
	private HashMap<String, char[]> seqlib;
	private ArrayList<String> targets;
	private HashMap<String, Integer> idToIndex = new HashMap<String, Integer>();
	private double[][] alignments;
	private Aligner gotoh;
	private BufferedWriter out = null;
	private BufferedReader tReader = null;
	private String targetsfolder;

	// family recognition test
	private int fam_recsamefam;
	private int fam_recsamesup;
	private int fam_recsamefold;
	private int fam_recdiffold;

	// supfam recognition
	private int supfam_recsamesup;
	private int supfam_recsamefold;
	private int supfam_recdiffold;

	// fold recognition
	private int fold_recsamefold;
	private int fold_recdiffold;

	public ValidateAlignments(Aligner gotoh, HashMap<String, char[]> sl,
			ArrayList<String> targets, String targetsfolder,
			BufferedWriter resultwriter) {
		this.gotoh = gotoh;
		this.seqlib = sl;
		this.targets = targets;
		alignments = new double[seqlib.size()][seqlib.size()];
		out = resultwriter;
		this.targetsfolder = targetsfolder;
		init();
	}

	private void init() {
		for (int i = 0; i < alignments[0].length; i++) {
			for (int j = 0; j < alignments.length; j++) {
				alignments[i][j] = Double.NEGATIVE_INFINITY;
			}
		}
	}

	public void benchmark() throws IOException {
		double score = Double.NEGATIVE_INFINITY;
		int besthit = 0;
		SequenceAlignment temp = null;

		int index = 0;
		int id1;
		int id2;

		String line;
		ArrayList<String[]> templates = new ArrayList<String[]>();

		for (String target : targets) {
			templates.removeAll(templates);
			double maxscore = Double.NEGATIVE_INFINITY;
			double sup_maxscore = Double.NEGATIVE_INFINITY;
			double fold_maxscore = Double.NEGATIVE_INFINITY;
			double outer_maxscore = Double.NEGATIVE_INFINITY;
			besthit = -1;
			String[] templevels;
			ArrayList<String> levels = new ArrayList<String>();		
			
			tReader = new BufferedReader(new FileReader(targetsfolder + "/"
					+ target + ".templates"));
			while ((line = tReader.readLine()) != null) {
				if (!line.startsWith("#")) {
					templates.add(line.split("\\s+"));
				} else {
					templevels = line.split("\\s+");
					for(String l : templevels){
						levels.add(l);
					}
				}
			}

			for (int i = 0; i < templates.size(); i++) {
				String template = templates.get(i)[0];
				String level = templates.get(i)[2];

				// give index to ID for storing the score of the alignment in a
				// matrix
				if (!idToIndex.containsKey(target)) {
					idToIndex.put(target, index);
					index++;
				}

				if (!idToIndex.containsKey(templates.get(i)[0])) {
					idToIndex.put(template, index);
					index++;
				}

				// make Alignment and check for besthit
				id1 = idToIndex.get(target);
				id2 = idToIndex.get(template);
				Sequence seq1 = new Sequence(target,seqlib.get(target));
				Sequence seq2 = new Sequence(template, seqlib.get(template));
				if (alignments[id1][id2] == Double.NEGATIVE_INFINITY) {
					temp = (SequenceAlignment) gotoh.align(seq1,seq2);
//					int minlength = Math.min(seq1.length(), seq2.length());
					score = (temp.getScore()/temp.countAlignedResidues());
					alignments[id1][id2] = score;
					alignments[id2][id1] = score;
				}		
				// else get score from alignment matrix
				else {
					score = alignments[id1][id2];
					
				}
//				out.append(target+"\t"+template+"\t"+((Math.round(score*100.0))/100.0)+"\n");
				if (score > maxscore) {
					maxscore = score;
					besthit = i;
				} 
				if (level.equals("supfam") && score > sup_maxscore) {
					sup_maxscore = score;
				} else if (level.equals("fold") && score > fold_maxscore) {
					fold_maxscore = score;
				} else if (level.equals("outer") && score > outer_maxscore) {
					outer_maxscore = score;
				}
			}

			// famrec test
			if(levels.contains("outer") && (levels.contains("fam") || levels.contains("supfam") || levels.contains("fold"))){
				if (templates.get(besthit)[2].equals("fam")) {
					fam_recsamefam++;
				} else if (templates.get(besthit)[2].equals("supfam")) {
					fam_recsamesup++;
				} else if (templates.get(besthit)[2].equals("fold")) {
					fam_recsamefold++;
				} else {
					fam_recdiffold++;
				}
			}

			// suprec test
			if((levels.contains("supfam")||levels.contains("fold")) && levels.contains("outer")){
				if (sup_maxscore > fold_maxscore) {
					if (sup_maxscore > outer_maxscore) {
						supfam_recsamesup++;
					} else {
						supfam_recdiffold++;
					}
				} else {
					if (fold_maxscore > outer_maxscore) {
						supfam_recsamefold++;
						
					} else {
						supfam_recdiffold++;
					}
				}
			}
			
			if(levels.contains("fold") && levels.contains("outer")){
				//foldrec test
				if(fold_maxscore > outer_maxscore){
					fold_recsamefold++;
				} else {
					fold_recdiffold++;
				}
			}
			System.out.println("found max");
		}
//		out.close();
	}

	// private void printResultList(CathScopEntry besthit, CathScopEntry query)
	// throws IOException {
	// out.append(query.getID()+"\t"+query.getCathClazz()+"."+query.getCathFold()+"."+query.getCathSupFam()+"."+query.getCathFam()+"\t");
	// out.append(query.getScopClazz()+"."+query.getScopFold()+"."+query.getScopSupFam()+"."+query.getScopFam()+"\t");
	// out.append(besthit.getID()+"\t"+besthit.getCathClazz()+"."+besthit.getCathFold()+"."+besthit.getCathSupFam()+"."+besthit.getCathFam()+"\t");
	// out.append(besthit.getScopClazz()+"."+besthit.getScopFold()+"."+besthit.getScopSupFam()+"."+besthit.getScopFam()+"\t"+alignments[idToIndex.get(query.getID())][idToIndex.get(besthit.getID())]+"\n");
	// }

	public void printResults() throws IOException {
		// fam recognition
		out.write(fam_recsamefam + "\n" + fam_recsamesup + "\n"
				+ fam_recsamefold + "\n" + fam_recdiffold + "\n");

		// supfam recognition
		out.write(supfam_recsamesup + "\n" + supfam_recsamefold + "\n"
				+ supfam_recdiffold + "\n");

		// fold recognition
		out.write(fold_recsamefold + "\n" + fold_recdiffold + "\n");

		out.close();
	}

	// public int[] getResults(){
	// int cath_foldrecognition =
	// cath_recsamefam+cath_recsamesup+cath_recsamefold;
	// int scop_foldrecognition =
	// scop_recsamefam+scop_recsamesup+scop_recsamefold;
	// return new int[]{cath_foldrecognition,scop_foldrecognition};
	// }

	public static void main(String[] args) throws IOException {
		double go = Double.parseDouble(args[0]);
		double ge = Double.parseDouble(args[1]);
		double hbWeight = Double.parseDouble(args[7]);
		double polWeight = Double.parseDouble(args[6]);
		double seqWeight = Double.parseDouble(args[5]);
		double ssWeight = Double.parseDouble(args[8]);
		HashMap<String, char[]> sec_seqlib = SeqLibrary
				.read(args[3]);
		HashMap<String,char[]> seqlib = SeqLibrary.read(args[2]);
		InitClass init = new InitClass();
		
		HashMap<String, Integer> t = new HashMap<String,Integer>();
		File folder = new File(args[10]);
		File[] files = folder.listFiles();
		for(File f : files){
			t.put(f.getName().split("\\.")[0],0);
		}
		BufferedReader in = new BufferedReader(new FileReader(args[4]));
		String line;
		ArrayList<String> targets = new ArrayList<String>();
		while((line = in.readLine())!= null){
			if(t.containsKey(line))	targets.add(line);
		}
		in.close();
		

		
//		for (int go = -5; go > -13; go--) {
//			for (int ge = -1; ge >= go ; ge--) {
				GlobalMusterLite muster = new GlobalMusterLite(go, ge, sec_seqlib,
						init.calcGotohInputMatrix(init.calcHydropathyScores()),
						init.calcGotohInputMatrix(init.calcPolarityScores()),
						SecStructScores.matrix, QuasarMatrix.DAYHOFF_MATRIX, hbWeight, polWeight,
						ssWeight, seqWeight);
				BufferedWriter resultWriter = new BufferedWriter(new FileWriter(args[9]));
				ValidateAlignments validate = new ValidateAlignments(muster, seqlib, targets, args[10], resultWriter);
				validate.benchmark();
				validate.printResults();
//			}
//		}
	}
}
