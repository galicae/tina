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
import bioinfo.proteins.PDBFileReader;
import bioinfo.superpos.TMMain;
import bioinfo.superpos.Transformation;

public class BMCathAndScop {

	private HashMap<String, char[]> seqlib;
	private HashMap<String, CathScopEntry> cathscopinfos;
	private HashMap<String, ArrayList<String>> pairs;

	private String pdbFolder;

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

	private Aligner gotoh;

	private BufferedWriter out = null;
	private BufferedWriter FamRecWr = null;
	private BufferedWriter scopFamMisWr = null;
	private BufferedWriter scopSupMisWr = null;
	private BufferedWriter scopFolMisWr = null;
	private BufferedWriter cathFamMisWr = null;
	private BufferedWriter cathSupMisWr = null;
	private BufferedWriter cathFolMisWr = null;


	public BMCathAndScop(String seqlibpath, String cathscopinfopath,
			String pairlistpath, String matrixpath, int go, int ge,
			String mode, String outpath, String pdbFolder, String famRecPath,
			String scopFamMisPath, String scopSupMisPath,
			String scopFolMisPath, String cathFamMisPath,
			String cathSupMisPath, String cathFolMisPath) {
		try {
			this.out = new BufferedWriter(new FileWriter(outpath));
		} catch (IOException e) {
			System.out.println("geht nischt");
		}
		this.seqlib = SeqLibrary.parse(seqlibpath);
		this.cathscopinfos = CathScopHash.read(cathscopinfopath);
		this.pairs = HashPairReader.readPairs(pairlistpath
				+ "/cathscop.inpairs", pairlistpath + "/cathscop.outpairs");
		double[][] scoringmatrix = QuasarMatrix.parseMatrix(matrixpath);
		this.pdbFolder = pdbFolder;

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

		if (mode.equals("global")) {
			gotoh = new GlobalSequenceGotoh(go, ge, scoringmatrix);
		} else if (mode.equals("local")) {
			gotoh = new LocalSequenceGotoh(go, ge, scoringmatrix);
		} else {
			gotoh = new FreeshiftSequenceGotoh(go, ge, scoringmatrix);
		}

		try {
			FamRecWr = new BufferedWriter(new FileWriter(famRecPath));
			scopFamMisWr = new BufferedWriter(new FileWriter(scopFamMisPath));
			scopSupMisWr = new BufferedWriter(new FileWriter(scopSupMisPath));
			scopFolMisWr = new BufferedWriter(new FileWriter(scopFolMisPath));
			cathFamMisWr = new BufferedWriter(new FileWriter(cathFamMisPath));
			cathSupMisWr = new BufferedWriter(new FileWriter(cathSupMisPath));
			cathFolMisWr = new BufferedWriter(new FileWriter(cathFolMisPath));
		} catch (Exception e) {
			e.printStackTrace();
		}

		benchmark();
	}

	public void benchmark() {
		PDBFileReader nikoReader = new PDBFileReader();
		nikoReader.setFolder(pdbFolder);
		TMMain tmMachine = new TMMain();

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
		
		int count = 0;
		for(Entry<String, ArrayList<String>> query_entry : pairs.entrySet()){
			count += query_entry.getValue().size();
		}
		System.out.println(count);
		for (Entry<String, ArrayList<String>> query_entry : pairs.entrySet()) {
			CathScopEntry besthit = null;
			maxscore = Double.NEGATIVE_INFINITY;
			samefammaxscore_cath = Double.NEGATIVE_INFINITY;
			samesupmaxscore_cath = Double.NEGATIVE_INFINITY;
			samefoldmaxscore_cath = Double.NEGATIVE_INFINITY;
			diffoldmaxscore_cath = Double.NEGATIVE_INFINITY;
			samefammaxscore_scop = Double.NEGATIVE_INFINITY;
			samesupmaxscore_scop = Double.NEGATIVE_INFINITY;
			samefoldmaxscore_scop = Double.NEGATIVE_INFINITY;
			diffoldmaxscore_scop = Double.NEGATIVE_INFINITY;
			SequenceAlignment famRec = null, cathFamMis = null, cathSupMis = null, cathFoldMis = null, cathDifFol = null; // temporary
																															// placeholders
																															// for
																															// max
																															// scoring
																															// alignments
																															// in
																															// corresponding
																															// benchmark
																															// strategies
			SequenceAlignment scopFamMis = null, scopSupMis = null, scopFoldMis = null, scopDifFol = null;
			SequenceAlignment temp;
			CathScopEntry query = cathscopinfos.get(query_entry.getKey());
			CathScopEntry template;
			
			for (String template_id : query_entry.getValue()) {
				temp = (SequenceAlignment) gotoh.align(
						new Sequence(query.getID(), seqlib.get(query.getID())),
						new Sequence(template_id, seqlib.get(template_id)));
				template = cathscopinfos.get(template_id);
				score = temp.getScore();
				if (score > maxscore) {
					maxscore = score;
					besthit = cathscopinfos.get(template_id);
					famRec = temp;
				}
				
				// for cath misclassification tests			
				if (query.getCathClazz() == template.getCathClazz()
						&& query.getCathFold() == template.getCathFold()
						&& query.getCathSupFam() == template.getCathSupFam()
						&& query.getCathFam() == template.getCathFam()) {
					if (score > samefammaxscore_cath) {
						samefammaxscore_cath = score;
						cathFamMis = temp;
					}
				}
				if (query.getCathClazz() == template.getCathClazz()
						&& query.getCathFold() == template.getCathFold()
						&& query.getCathSupFam() == template.getCathSupFam() 
						&& query.getCathFam() != template.getCathFam()) {
					if (score > samesupmaxscore_cath) {
						samesupmaxscore_cath = score;
						cathSupMis = temp;
					}
				}
				if (query.getCathClazz() == template.getCathClazz()
						&& query.getCathFold() == template.getCathFold()
						&& query.getCathSupFam() != template.getCathSupFam()) {
					if (score > samefoldmaxscore_cath) {
						samefoldmaxscore_cath = score;
						cathFoldMis = temp;
					}
				}
				if (query.getCathClazz() != template.getCathClazz()
						|| query.getCathFold() != template.getCathFold()) {
					if (score > diffoldmaxscore_cath) {
						diffoldmaxscore_cath = score;
						cathDifFol = temp;
					}
				}

				// for scop misclassification tests
				if (query.getScopClazz() == template.getScopClazz()
						&& query.getScopFold() == template.getScopFold()
						&& query.getScopSupFam() == template.getScopSupFam()
						&& query.getScopFam() == template.getScopFam()) {
					if (score > samefammaxscore_scop) {
						samefammaxscore_scop = score;
						scopFamMis = temp;
					}
				}
				if (query.getScopClazz() == template.getScopClazz()
						&& query.getScopFold() == template.getScopFold()
						&& query.getScopSupFam() == template.getScopSupFam()
						&& query.getScopFam() != template.getScopFam()) {
					if (score > samesupmaxscore_scop) {
						samesupmaxscore_scop = score;
						scopSupMis = temp;
					}
				}
				if (query.getScopClazz() == template.getScopClazz()
						&& query.getScopFold() == template.getScopFold()
						&& query.getScopSupFam() != template.getScopSupFam()) {
					if (score > samefoldmaxscore_scop) {
						samefoldmaxscore_scop = score;
						scopFoldMis = temp;
					}
				}
				if (query.getScopClazz() != template.getScopClazz()
						|| query.getScopFold() != template.getScopFold()) {
					if (score > diffoldmaxscore_scop) {
						diffoldmaxscore_scop = score;
						scopDifFol = temp;
					}
				}
			}

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
			try {
				Transformation tr = tmMachine.calculateTransformation(famRec,
						nikoReader.readFromFolderById(famRec.getComponent(0)
								.getID()), nikoReader.readFromFolderById(famRec
								.getComponent(1).getID()));
				FamRecWr.write(tr.getTmscore() + " " + tr.getGdt() + "\n");
			} catch (Exception e) {
				e.printStackTrace();
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

			// cath misclassification test
			if (samefammaxscore_cath >= diffoldmaxscore_cath && diffoldmaxscore_cath != Double.NEGATIVE_INFINITY) {
				cath_missamefam++;
				try {
					Transformation tr = tmMachine.calculateTransformation(
							cathFamMis, nikoReader
									.readFromFolderById(cathFamMis
											.getComponent(0).getID()),
							nikoReader.readFromFolderById(cathFamMis
									.getComponent(1).getID()));
					cathFamMisWr.write(tr.getTmscore() + " " + tr.getGdt()
							+ "\n");
				} catch (Exception e) {
					e.printStackTrace();
				}
			} else if(samefammaxscore_cath < diffoldmaxscore_cath && samefammaxscore_cath != Double.NEGATIVE_INFINITY) {
				cath_misdiffoldfam++;
				try {
					Transformation tr = tmMachine.calculateTransformation(
							cathDifFol, nikoReader
									.readFromFolderById(cathDifFol
											.getComponent(0).getID()),
							nikoReader.readFromFolderById(cathDifFol
									.getComponent(1).getID()));
					cathFamMisWr
							.write(tr.getTmscore() + " " + tr.getGdt() + "\n");
				} catch (Exception e) {
					e.printStackTrace();
				}
			}
			if (samesupmaxscore_cath >= diffoldmaxscore_cath && diffoldmaxscore_cath != Double.NEGATIVE_INFINITY) {
				cath_missamesup++;
				try {
					Transformation tr = tmMachine.calculateTransformation(
							cathSupMis, nikoReader
									.readFromFolderById(cathSupMis
											.getComponent(0).getID()),
							nikoReader.readFromFolderById(cathSupMis
									.getComponent(1).getID()));
					cathSupMisWr.write(tr.getTmscore() + " " + tr.getGdt()
							+ "\n");
				} catch (Exception e) {
					e.printStackTrace();
				}
			} else if(samesupmaxscore_cath < diffoldmaxscore_cath && samesupmaxscore_cath != Double.NEGATIVE_INFINITY){
				cath_misdiffoldsup++;
				try {
					Transformation tr = tmMachine.calculateTransformation(
							cathDifFol, nikoReader
									.readFromFolderById(cathDifFol
											.getComponent(0).getID()),
							nikoReader.readFromFolderById(cathDifFol
									.getComponent(1).getID()));
					cathSupMisWr
							.write(tr.getTmscore() + " " + tr.getGdt() + "\n");
				} catch (Exception e) {
					e.printStackTrace();
				}
			}
			if (samefoldmaxscore_cath >= diffoldmaxscore_cath && diffoldmaxscore_cath != Double.NEGATIVE_INFINITY) {
				cath_missamefold++;
				try {
					Transformation tr = tmMachine.calculateTransformation(
							cathFoldMis, nikoReader
									.readFromFolderById(cathFoldMis
											.getComponent(0).getID()),
							nikoReader.readFromFolderById(cathFoldMis
									.getComponent(1).getID()));
					cathFolMisWr.write(tr.getTmscore() + " " + tr.getGdt()
							+ "\n");
				} catch (Exception e) {
					e.printStackTrace();
				}
			} else if(samefoldmaxscore_cath < diffoldmaxscore_cath && samefoldmaxscore_cath != Double.NEGATIVE_INFINITY) {
				cath_misdiffoldfold++;
				try {
					Transformation tr = tmMachine.calculateTransformation(
							cathDifFol, nikoReader
									.readFromFolderById(cathDifFol
											.getComponent(0).getID()),
							nikoReader.readFromFolderById(cathDifFol
									.getComponent(1).getID()));
					cathFolMisWr
							.write(tr.getTmscore() + " " + tr.getGdt() + "\n");
				} catch (Exception e) {
					e.printStackTrace();
				}
			}

			// scop misclassification test
			if (samefammaxscore_scop >= diffoldmaxscore_scop && diffoldmaxscore_scop != Double.NEGATIVE_INFINITY) {
				scop_missamefam++;
				try {
					Transformation tr = tmMachine.calculateTransformation(
							scopFamMis, nikoReader
									.readFromFolderById(scopFamMis
											.getComponent(0).getID()),
							nikoReader.readFromFolderById(scopFamMis
									.getComponent(1).getID()));
					scopFamMisWr.write(tr.getTmscore() + " " + tr.getGdt()
							+ "\n");
				} catch (Exception e) {
					e.printStackTrace();
				}
			} else if(samefammaxscore_scop < diffoldmaxscore_scop && samefammaxscore_scop != Double.NEGATIVE_INFINITY) {
				scop_misdiffoldfam++;
				try {
					Transformation tr = tmMachine.calculateTransformation(
							scopDifFol, nikoReader
									.readFromFolderById(scopDifFol
											.getComponent(0).getID()),
							nikoReader.readFromFolderById(scopDifFol
									.getComponent(1).getID()));
					scopFamMisWr
							.write(tr.getTmscore() + " " + tr.getGdt() + "\n");
				} catch (Exception e) {
					e.printStackTrace();
				}
			}
			if (samesupmaxscore_scop >= diffoldmaxscore_scop && diffoldmaxscore_scop != Double.NEGATIVE_INFINITY) {
				scop_missamesup++;
				try {
					Transformation tr = tmMachine.calculateTransformation(
							scopSupMis, nikoReader
									.readFromFolderById(scopSupMis
											.getComponent(0).getID()),
							nikoReader.readFromFolderById(scopSupMis
									.getComponent(1).getID()));
					scopSupMisWr.write(tr.getTmscore() + " " + tr.getGdt()
							+ "\n");
				} catch (Exception e) {
					e.printStackTrace();
				}
			} else if(samesupmaxscore_scop < diffoldmaxscore_scop && samesupmaxscore_scop != Double.NEGATIVE_INFINITY) {
				scop_misdiffoldsup++;
				try {
					Transformation tr = tmMachine.calculateTransformation(
							scopDifFol, nikoReader
									.readFromFolderById(scopDifFol
											.getComponent(0).getID()),
							nikoReader.readFromFolderById(scopDifFol
									.getComponent(1).getID()));
					scopSupMisWr
							.write(tr.getTmscore() + " " + tr.getGdt() + "\n");
				} catch (Exception e) {
					e.printStackTrace();
				}
			}
			if (samefoldmaxscore_scop >= diffoldmaxscore_scop && diffoldmaxscore_scop != Double.NEGATIVE_INFINITY) {
				scop_missamefold++;
				try {
					Transformation tr = tmMachine.calculateTransformation(
							scopFoldMis, nikoReader
									.readFromFolderById(scopFoldMis
											.getComponent(0).getID()),
							nikoReader.readFromFolderById(scopFoldMis
									.getComponent(1).getID()));
					scopFolMisWr.write(tr.getTmscore() + " " + tr.getGdt()
							+ "\n");
				} catch (Exception e) {
					e.printStackTrace();
				}
			} else if(samefoldmaxscore_scop < diffoldmaxscore_scop && samefoldmaxscore_scop != Double.NEGATIVE_INFINITY) {
				scop_misdiffoldfold++;
				try {
					Transformation tr = tmMachine.calculateTransformation(
							scopDifFol, nikoReader
									.readFromFolderById(scopDifFol
											.getComponent(0).getID()),
							nikoReader.readFromFolderById(scopDifFol
									.getComponent(1).getID()));
					scopFolMisWr
							.write(tr.getTmscore() + " " + tr.getGdt() + "\n");
				} catch (Exception e) {
					e.printStackTrace();
				}
			}
			System.out.println("found max for id: "+query_entry.getKey());
		}

		try {
			printResults();
			FamRecWr.close();
			scopFamMisWr.close();
			scopSupMisWr.close();
			scopFolMisWr.close();
			cathFamMisWr.close();
			cathSupMisWr.close();
			cathFolMisWr.close();
			out.close();
		} catch (IOException e) {
			System.out.println("Cannot write results");
		}
		System.out.println("ready");
	}

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
	}

	public static void main(String[] args) {
		BMCathAndScop benchmark = new BMCathAndScop(args[0], args[1], args[2],
				args[3], Integer.parseInt(args[4]), Integer.parseInt(args[5]),
				args[6], args[7], args[8], args[9], args[10], args[11],
				args[12], args[13], args[14], args[15]);
	}
}
