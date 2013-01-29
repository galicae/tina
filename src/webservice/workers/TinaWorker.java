/******************************************************************************
 * webservice.workers.TinaWorker                                              *
 *                                                                            *
 * This file is best read at line width 80 and tab width 4.                   *
 ******************************************************************************/
package webservice.workers;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.LinkedList;
import bioinfo.Sequence;
import bioinfo.alignment.SequenceAlignment;
import bioinfo.alignment.gotoh.FreeshiftSequenceGotoh;
import bioinfo.pdb.PDBFile;
import bioinfo.proteins.PDBEntry;
import bioinfo.proteins.PDBFileReader;
import bioinfo.proteins.structure.SimpleCoordMapper;

/**
 * TinaWorker is the worker that works on an TINA job.
 * @author huberste
 * @date December 09, 2012
 * @version alpha
 */
public class TinaWorker extends Worker {
	
	private final static String PDB_FILE_PATH = "/home/h/huberste/gobi/webserver/pdb/";
	
	private final static String SEQLIB_FILE = "/home/h/huberste/gobi/webserver/database/domains.seqlib";
	private final static double GO = -10.0;
	private final static double GE = -2.0;
	private final static double[][] MATRIX = bioinfo.alignment.matrices.QuasarMatrix.DAYHOFF_MATRIX;
	
	private final static int TOP_NUM = 5;
	
	private final static Sequence BAD = new Sequence("id", "A");
	private final static SequenceAlignment WORST = new SequenceAlignment(
					BAD, BAD, "A","A", Double.NEGATIVE_INFINITY);
	
	private Sequence sequence;
	private LinkedList<Sequence> sequencedb;
	
	String result;
	
	/**
	 * 
	 * @param jobfile
	 */
	public TinaWorker(String jobfile) {
		super(jobfile);
	}
	
	/**
	 * This will make a prediction of the structure.
	 * The plan is:
	 * 1) align with freeshiftGotoh against all Sequences in $SEQLIB_FILE.
	 * 2) take (five?) best alignment(s). for each:
	 * 3) Thread against these alignments.
	 */
	public void work() {
		// Read Sequence Database
		getSequenceDB();
		readFile();
		
		// DONE 1) Align against all sequences in $SEQLIB_FILE
		// && DONE 2) take TOP_NUM best Alignments
		SequenceAlignment[] topAlis = new SequenceAlignment[TOP_NUM];
		for (int i = 0; i < TOP_NUM; i++) {
			topAlis[i]=WORST;
		}
		// smallest score in top TOP_NUM
		double leastScore = Double.NEGATIVE_INFINITY;
		// smallest TOP_NUM
		int leastScoring = 0;
		FreeshiftSequenceGotoh gotoh = new FreeshiftSequenceGotoh(GO, GE, MATRIX);
		for (Sequence seq : sequencedb) {
			// DONE debugging: what are the two sequences???
//			System.out.println("debugging: work() l. 83: seq = "+seq.toStringVerbose());
//			System.out.println("debugging: work() l. 83: sequence= "+sequence.toStringVerbose());
			SequenceAlignment temp = gotoh.align(sequence, seq);
			// DONE debugging. Am I here?
//			System.out.println("debugging: work() l. 83: temp = "+temp.toStringVerbose());
			if (temp.getScore() > leastScore) {
				topAlis[leastScoring] = temp;
				// calculate next smallest
				leastScore = topAlis[0].getScore();
				leastScoring = 0;
				for (int i = 0; i < TOP_NUM; i++) {
					if (topAlis[i].getScore() < leastScore) {
						leastScore = topAlis[i].getScore();
						leastScoring = i;
					}
				}
			}
		}
		
		BufferedReader from = null;
		BufferedWriter to = null;
		
		String line = null; 
		try {
			from = new BufferedReader(new FileReader(JOB_FILE));
			to = new BufferedWriter(new FileWriter(DONE_FILE));
			while((line = from.readLine()) != null) {
				to.write(line+"\n");
			}
			to.write("RESULT=\n");
//			to.write(result);
			
			to.write("TOP_ALIGNMENTS=\n");
			for (int i = 0; i < TOP_NUM; i++) {
				to.write(topAlis[i].toStringVerbose()+"\n");
			}
			
			// 3) Thread Alignments
			// get PDBFile(s) of top5
			PDBFileReader pdbReader = new PDBFileReader(PDB_FILE_PATH);
			PDBEntry[] pdbentries = new PDBEntry[TOP_NUM];
			for (int i = 0; i < TOP_NUM; i++) {
				String idi = topAlis[i].getComponent(1).getID().toUpperCase();
				// Get PDB Files from PDB Server
				PDBFile.getFile(PDB_FILE_PATH, idi.substring(0, 4));
				// read Entrys from PDBFiles
				pdbentries[i] = pdbReader.readFromFolderById(idi);
			}
			
			// Thread topAlis[i] with pdbentries[i]
			
			PDBEntry[] results = new PDBEntry[TOP_NUM];
			
			for (int i = 0; i < TOP_NUM; i++) {
				results[i] = SimpleCoordMapper.map(topAlis[i], pdbentries[i]);
			}
			
			to.write("PDBs=\n");
			// write results in result String
			for (int i = 0; i < TOP_NUM; i++) {
				
				// TODO write better ToString function!
//				to.write(results[i].toString()+"\n");
				to.write(results[i].getAtomSectionAsString());
			}
			
			} catch (IOException e) {
				System.err.println("Error while trying to copy "+JOB_FILE+" to "+DONE_FILE+".");
				e.printStackTrace();
			} finally {
				try {
					if (from != null) from.close();
					if (to != null) to.close();
				} catch (IOException e) {
					System.err.println("Error while trying close FileStreams");
					e.printStackTrace();
				}
			}
		// rm $JOBS_DIR/02_working/$JOB_ID.job
		new File(JOB_FILE).delete();
		
		// save results in $JOB_DIR/03_done/$JOB_ID.job
		// and rm $JOB_DIR/02_working/$JOB_ID.job
//		writeResult();
	}

	/**
	 * Loads the Sequences in the given $SEQLIB_FILE into LinkedList $sequencedb
	 */
	private void getSequenceDB() {
		sequencedb = new LinkedList<Sequence>();

		BufferedReader reader = null;
		try {
			reader = new BufferedReader(new FileReader(SEQLIB_FILE));
			String line = null;
			while ((line = reader.readLine()) != null) {
				// TODO optimization
				String[] temp= line.split(":");
				sequencedb.add(new Sequence(temp[0], temp[1].toCharArray()));
			}
		} catch (IOException x) {
			System.err.format("IOException: %s%n", x);
		} finally {
			try {
				reader.close();
			} catch (IOException x) {
				System.err.format("IOException: %s%n", x);
			}
		}
	}
	
	/**
	 * reads the jobfile and all important data from it.
	 */
	@Override
	protected void readFile() {
		BufferedReader from = null;
		
		String line = null; 
		try {
			from = new BufferedReader(new FileReader(JOB_FILE));
			while((line=from.readLine())!= null) {
				if (line.startsWith("SEQUENCE=")) {
					String[] temp = line.substring(9).split(":");
					// DONE stupid user: 
					if (temp.length < 2) { // Sequence not given in Format "id:sequence"
						result = "INPUT ERROR: SEQUENCE ONE WAS NOT GIVEN IN FORMAT \"ID:SEQUENCE\"";
					}
					sequence = new Sequence(temp[0].trim(), temp[1].trim());
					// DONE debugging
					//System.out.println("debugging: sequence = "+sequence.toStringVerbose());
				}
			}
		} catch (IOException e) {
			System.err.println("Error while trying to read "+JOB_FILE+".");
			e.printStackTrace();
		} finally {
			try {
				from.close();
			} catch (IOException e) {
				System.err.println("Error while trying close "+JOB_FILE+".");
				e.printStackTrace();
			}
		}
		
	}

	/**
	 * Writes the result file.
	 * Also deletes the working file.
	 */
	@Override
	protected void writeResult() {
		BufferedReader from = null;
		BufferedWriter to = null;
		
		String line = null; 
		try {
			from = new BufferedReader(new FileReader(JOB_FILE));
			to = new BufferedWriter(new FileWriter(DONE_FILE));
			while((line = from.readLine()) != null) {
				to.write(line+"\n");
			}
			to.write("RESULT=\n");
			to.write(result);
		} catch (IOException e) {
			System.err.println("Error while trying to copy "+JOB_FILE+" to "+DONE_FILE+".");
			e.printStackTrace();
		} finally {
			try {
				if (from != null) from.close();
				if (to != null) to.close();
			} catch (IOException e) {
				System.err.println("Error while trying close FileStreams");
				e.printStackTrace();
			}
		}
		// rm $JOBS_DIR/02_working/$JOB_ID.job
		new File(JOB_FILE).delete();
	}
	
}