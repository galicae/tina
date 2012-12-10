/**
 * 
 */
package webservice.workers;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;

import bioinfo.Sequence;
import bioinfo.alignment.SequenceAlignment;
import bioinfo.alignment.gotoh.FreeshiftSequenceGotoh;
import bioinfo.alignment.matrices.QuasarMatrix;

/**
 * @author huberste
 *
 */
public class GotohWorker extends Worker{

	private Sequence sequenceOne;
	private Sequence sequenceTwo;
	boolean shortMatrix;
	double[][] matrix;
	
	SequenceAlignment aliresult;
	String result;
	
	public GotohWorker(String jobFile) {
		super(jobFile);
	}

	@Override
	public void work() {
		// read WORKING job file
		readFile();
		
		// DONE debugging: Check if matrix was given correctly.
//			String matrixString = "";
//			for (int x = 0; x < matrix.length; x++) {
//				matrixString += "(";
//				for (int y = 0; y < matrix[x].length; y++) {
//					matrixString += matrix[x][y] + ", ";
//				}
//				matrixString += ")";
//			}
//			System.out.println("Matrix:\n"+matrixString);
		// end debugging
		
		FreeshiftSequenceGotoh gotoh = new FreeshiftSequenceGotoh(-10.0, -2.0, matrix);
		// TODO debugging: What are sequenceone and sequence two?
			System.out.println("SequenceOne = "+sequenceOne.toStringVerbose());
			System.out.println("SequenceTwo = "+sequenceTwo.toStringVerbose());
		// end debugging
		aliresult = gotoh.align(sequenceOne, sequenceTwo);
		
		// Write DONE job file
		writeResult();
	}

	@Override
	protected void readFile() {
		BufferedReader from = null;
		
		String line = null; 
		try {
			from = new BufferedReader(new FileReader(JOB_FILE));
			while((line=from.readLine())!= null) {
				if (line.startsWith("SEQUENCE_ONE=")) {
					String[] temp = line.substring(13).split(":");
					// TODO stupid user: 
					if (temp.length < 2) { // Sequence not given in Format "id:sequence"
						result = "INPUT ERROR: SEQUENCE ONE WAS NOT GIVEN IN FORMAT \"ID:SEQUENCE\"";
					}
					sequenceOne = new Sequence(temp[0], temp[1]);
					// TODO debugging
					System.out.println("debugging: sequenceOne = "+sequenceOne.toStringVerbose());
				} else if (line.startsWith("SEQUENCE_TWO=")) {
					String[] temp = line.substring(13).split(":");
					if (temp.length < 2) { // Sequence not given in Format "id:sequence"
						result = "INPUT ERROR: SEQUENCE TWO WAS NOT GIVEN IN FORMAT \"ID:SEQUENCE\"";
					}
					sequenceTwo = new Sequence(temp[0], temp[1]);
					// TODO debugging
					System.out.println("debugging: sequenceTwo = "+sequenceTwo.toStringVerbose());
				} else if (line.startsWith("SHORT_MATRIX=")) {
					String temp = line.substring(13);
					shortMatrix = Boolean.parseBoolean(temp);
					// TODO debugging
					System.out.println("debugging: shortMatrix = "+String.valueOf(shortMatrix));
				} else if (line.startsWith("MATRIX=")) {
					String temp = "";
					while((line=from.readLine())!= null) {
						temp += line+"\n";
					}
					if (temp.startsWith("dayhoff")) {
						matrix = QuasarMatrix.DAYHOFF_MATRIX;
					} else {
						// TODO other matrices
					}
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
			to.write(aliresult.toStringVerbose());
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
		new File(JOB_FILE).delete();
	}

}
