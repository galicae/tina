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
	
	public GotohWorker(String jobFile) {
		super(jobFile);
	}

	@Override
	public void work() {
		readFile();
		FreeshiftSequenceGotoh gotoh = new FreeshiftSequenceGotoh(-10.0, -2.0, matrix);
		SequenceAlignment alignment = gotoh.align(sequenceOne, sequenceTwo);
		
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
			to.write(alignment.toStringVerbose());
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

	@Override
	protected void readFile() {
		BufferedReader from = null;
		
		String line = null; 
		try {
			from = new BufferedReader(new FileReader(JOB_FILE));
			while((line=from.readLine())!= null) {
				if (line.startsWith("SEQUENCE_ONE=")) {
					String[] temp = line.substring(13).split(":");
					sequenceOne = new Sequence(temp[0], temp[1]);
				} else if (line.startsWith("SEQUENCE_TWO=")) {
					String[] temp = line.substring(13).split(":");
					sequenceTwo = new Sequence(temp[0], temp[1]);
				} else if (line.startsWith("SHORT_MATRIX=")) {
					String temp = line.substring(13);
					shortMatrix = Boolean.parseBoolean(temp);
				} else if (line.startsWith("MATRIX=")) {
					String temp = line.substring(7);
					shortMatrix = Boolean.parseBoolean(temp);
					while((line=from.readLine())!= null) {
						temp += "\n"+line;
					}
					if (matrix.equals("dayhoff")) {
						matrix = QuasarMatrix.DAYHOFF_MATRIX;
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
		// TODO Auto-generated method stub
		
	}

}
