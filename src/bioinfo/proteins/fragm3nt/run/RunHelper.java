package bioinfo.proteins.fragm3nt.run;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.ByteArrayOutputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.util.LinkedList;

import org.apache.commons.exec.CommandLine;
import org.apache.commons.exec.DefaultExecutor;
import org.apache.commons.exec.PumpStreamHandler;

import bioinfo.Sequence;
import bioinfo.alignment.SequenceAlignment;
import bioinfo.alignment.gotoh.FreeshiftSequenceGotoh;
import bioinfo.alignment.gotoh.GlobalSequenceGotoh;
import bioinfo.alignment.gotoh.LocalSequenceGotoh;
import bioinfo.alignment.matrices.QuasarMatrix;
import bioinfo.proteins.PDBEntry;
import bioinfo.proteins.PDBFileReader;
import bioinfo.proteins.fragm3nt.ProteinFragment;
import bioinfo.proteins.fragm3nt.assembler.AlignmentAssembler;
import bioinfo.superpos.TMMain;
import bioinfo.superpos.Transformation;

/**
 * this class is meant to be a static import in most of the classes in the run
 * package. It contains the functions most commonly used (creating/evaluating
 * Gotoh, doMagic(), fasta format, etc)
 * 
 * @author galicae
 * 
 */
public class RunHelper {

	static int fragLength = 8;
	static String desktop = "/home/galicae/Desktop/STRUCTURES/";
	static BufferedWriter resultWriter;

	/**
	 * this function calculates the gotoh part of things: for two PDB ids it
	 * creates the alignment with the best TM score and returns the TM score.
	 * 
	 * @param seq1ID
	 *            the id of the first sequence
	 * @param seq2ID
	 *            the id of the second sequence
	 * @return the optimal TM score of the superpositions of the proteins,
	 *         chosen from the best alignment (local, global, freeshift)
	 * @throws Exception
	 */
	public static double gotohPart(String seq1ID, String seq2ID)
			throws Exception {

		FreeshiftSequenceGotoh f = new FreeshiftSequenceGotoh(-9, -3,
				QuasarMatrix.DAYHOFF_MATRIX);
		GlobalSequenceGotoh g = new GlobalSequenceGotoh(-12, -1,
				QuasarMatrix.DAYHOFF_MATRIX);
		LocalSequenceGotoh l = new LocalSequenceGotoh(-12, -1,
				QuasarMatrix.DAYHOFF_MATRIX);

		Sequence seq1 = retrieveSeq(seq1ID);
		Sequence seq2 = retrieveSeq(seq2ID);

		SequenceAlignment fAli = f.align(seq1, seq2);
		SequenceAlignment gAli = g.align(seq1, seq2);
		SequenceAlignment lAli = l.align(seq1, seq2);

		BufferedWriter lw = new BufferedWriter(new FileWriter(
				"./fastaFiles/loc" + seq1.getID() + "_" + seq2.getID()));
		BufferedWriter fw = new BufferedWriter(new FileWriter(
				"./fastaFiles/fre" + seq1.getID() + "_" + seq2.getID()));
		BufferedWriter gw = new BufferedWriter(new FileWriter(
				"./fastaFiles/glo" + seq1.getID() + "_" + seq2.getID()));
		lw.write(toFastaFormat(lAli));
		fw.write(toFastaFormat(fAli));
		gw.write(toFastaFormat(gAli));
		lw.close();
		fw.close();
		gw.close();

		double fScore = 0;
		double gScore = 0;
		double lScore = 0;

		String fcall = "";
		String gcall = "";
		String lcall = "";

		String fOut = "";
		String gOut = "";
		String lOut = "";

		lcall = ("./tools/TMalign ");
		lcall += (desktop + seq1ID + ".pdb ");
		lcall += (desktop + seq2ID + ".pdb ");
		lcall += "-I ./fastaFiles/loc" + seq1ID + "_" + seq2ID;
		lOut = execToString(lcall);
		lScore = findMeTmScore(lOut);

		gcall = ("./tools/TMalign ");
		gcall += (desktop + seq1ID + ".pdb ");
		gcall += (desktop + seq2ID + ".pdb ");
		gcall += "-I ./fastaFiles/glo" + seq1ID + "_" + seq2ID;
		gOut = execToString(gcall);
		gScore = findMeTmScore(gOut);

		fcall = ("./tools/TMalign ");
		fcall += (desktop + seq1ID + ".pdb ");
		fcall += (desktop + seq2ID + ".pdb ");
		fcall += "-I ./fastaFiles/fre" + seq1ID + "_" + seq2ID;
		fOut = execToString(fcall);
		fScore = findMeTmScore(fOut);

		return Math.max(fScore, Math.max(gScore, lScore));
	}

	/**
	 * this function is the equivalent of the gotohPart for the
	 * AlignmentAssembler. It takes a list of sequences as input, performs the
	 * pairwise alignments from the sequences and feeds them to the
	 * AlignmentAssembler. The model produced from the method is checked against
	 * the native structure and the TM score of the optimal superposition is
	 * returned
	 * 
	 * @param id
	 *            the list of ids of the sequences deemed to be good templates
	 *            for the query protein
	 * @return the TM score of the superposition of the prediction of
	 *         AlignmentAssembler based on the alignments with the input
	 *         sequences
	 * @throws Exception
	 */
	public static double doMagic(LinkedList<String> id) throws Exception {
		// find all sequences in the id list and load them
		LinkedList<Sequence> seqs = loadSequences(id);

		// find all corresponding structures (ProteinFragment)
		LinkedList<PDBEntry> structures = loadPDBs(id);

		// now do predicting magic
		int extent = 5;
		AlignmentAssembler ass = new AlignmentAssembler(fragLength);
		ProteinFragment pred = ass.predStrucFromAl(seqs, extent, desktop);
		PDBEntry prediction = pred.toPDB();

		TMMain main = new TMMain();
		Sequence s = seqs.getFirst();
		char[] row = s.getSequence();
		SequenceAlignment alignment = new SequenceAlignment(s, s, row, row, 5);
		Transformation t = main.calculateTransformation(alignment,
				structures.get(0), prediction);
		return t.getTmscore();
	}

	public static LinkedList<PDBEntry> loadPDBs(LinkedList<String> ids) {
		LinkedList<PDBEntry> pdbs = new LinkedList<PDBEntry>();
		PDBFileReader reader = new PDBFileReader(desktop);
		for (int i = 0; i < ids.size(); i++) {
			pdbs.add(reader.readFromFolderById(ids.get(i)));
		}
		return pdbs;
	}

	/**
	 * this function takes a list of ids as input and returns the corresponding
	 * Sequence objects as a list. A cleaner implementation would have the
	 * directory as a parameter.
	 * 
	 * @param ids
	 *            the list of ids of sequences to return
	 * @return a list of the Sequence objects corresponding to the ids
	 */
	public static LinkedList<Sequence> loadSequences(LinkedList<String> ids) {
		LinkedList<Sequence> result = new LinkedList<Sequence>();
		try {
			BufferedReader r = new BufferedReader(new FileReader(
					"domains.seqlib"));
			String line = "";
			int count = 0;
			String stringSeq = "";
			while (count < ids.size() && (line = r.readLine()) != null) {
				for (int i = 0; i < ids.size(); i++) {
					if (line.startsWith(ids.get(i))) {
						stringSeq = line.split(":")[1];
						result.add(new Sequence(ids.get(i), stringSeq));
						count++;
						continue;
					}
				}
			}
			r.close();
		} catch (Exception e) {
			e.printStackTrace();
		}
		return result;
	}

	/**
	 * crude method of parsing the TMalign output tunneled here as a String
	 * 
	 * @param tm
	 *            the output of a TMalign run
	 * @return the TM score calculated by TMalign
	 */
	public static double findMeTmScore(String tm) {
		String[] r = tm.split("TM-score= ");
		String score = r[r.length - 1].substring(0, 7);
		return Double.parseDouble(score);
	}

	/**
	 * executes a command using the Apache exec library and returns the console
	 * output of the process as a String
	 * 
	 * @param command
	 *            the command line method call, including parameters and all
	 *            (really, command line copy-paste)
	 * @return the standard output of the process called
	 * @throws Exception
	 */
	public static String execToString(String command) throws Exception {
		ByteArrayOutputStream outputStream = new ByteArrayOutputStream();
		CommandLine commandline = CommandLine.parse(command);
		DefaultExecutor exec = new DefaultExecutor();
		PumpStreamHandler streamHandler = new PumpStreamHandler(outputStream);
		exec.setStreamHandler(streamHandler);
		exec.execute(commandline);
		return (outputStream.toString());
	}

	/**
	 * method to retrieve a single sequence out of domains.seqlib (again, a
	 * clean version would include the directory as a parameter, but it was late
	 * that night)
	 * 
	 * @param id
	 *            the PDB id of the sequence. Should be in domains.seqlib.
	 * @return the sequence object corresponding to the id
	 * @throws Exception
	 */
	public static Sequence retrieveSeq(String id) throws Exception {
		BufferedReader r = new BufferedReader(new FileReader("domains.seqlib"));
		String line = "";
		while ((line = r.readLine()) != null) {
			if (line.startsWith(id)) {
				String[] s = line.split(":");
				Sequence sq = new Sequence(s[0], s[1]);
				r.close();
				return sq;
			}
		}
		r.close();
		return null;
	}

	/**
	 * writes an alignment in FASTA format (that is, start with ">" and insert
	 * linebreaks every 80 characters)
	 * 
	 * @param ali
	 *            the alignment to turn to FASTA format
	 * @return
	 */
	public static String toFastaFormat(SequenceAlignment ali) {
		StringBuilder sb = new StringBuilder();
		sb.append(">" + ali.getComponent(0).getID() + "\n");
		for (int i = 0; i < ali.getRowAsString(0).length(); i++) {
			sb.append(ali.getRow(0)[i]);
			if (i % 80 == 0 && i > 0)
				sb.append("\n");
		}
		sb.append("\n");
		sb.append(">" + ali.getComponent(1).getID() + "\n");
		for (int i = 0; i < ali.getRowAsString(1).length(); i++) {
			sb.append(ali.getRow(1)[i]);
			if (i % 80 == 0 && i > 0)
				sb.append("\n");
		}
		return sb.toString();
	}
}
