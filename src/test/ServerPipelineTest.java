/**
 * 
 */
package test;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;

/**
 * @author huberste
 *
 */
public class ServerPipelineTest {

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		
//		final String PDB_FILE_PATH = "/home/h/huberste/gobi/webserver/pdb/";
		
		final String WEBAPP_HOME = "/home/h/huberste/gobi/webserver/webapp";
		final String JOBS_DIR = "/home/h/huberste/gobi/webserver/jobs";
		
		String[] arg = new String[4];
		arg[0] = JOBS_DIR;
		arg[1] = WEBAPP_HOME;
		arg[2] = "tina";
//		arg[3] = "1a2xB00:EEKRNRAITARRQHLKSVMLQIAATELEKE";
//		arg[3] = "1tim:ASDFASDFASDFASDFASDF";
//		arg[3] = "1bla:ARGHHARLKSIGDAS";
		arg[3] = "1v7uA00:LPAHGCRHVAIIMDGNGRWAKKQGKIRAFGHKAGAKSVRRAVSFAANNGIEALTLYAFSSENWNRPAQEVSALMELFVWALDSEVKSLHRHNVRLRIIGDTSRFNSRLQERIRKSEALTAGNTGLTLNIAANYGGRWDIVQGVRQLAEKVQQGNLQPDQIDEEMLNQHVCMHELAPVDLVIRTGGEHRISNFLLWQIAYAELYFTDVLWPDFDEQDFEGALNAFAN";
		System.out.println("Job created with ID "+startJob(arg)+".");
		
		// TODO debugging
			// Comment: Debugging the debugging file... debuggingception.
//		PDBFileReader pdbReader = new PDBFileReader(PDB_FILE_PATH);
//		String idi = "1a2xa00";
//		PDBEntry temp = pdbReader.readFromFolderById(idi.toUpperCase());
//		System.out.println(temp);
		
	}
	
	/**
	 * starts a new Job in this process. For testing purpose only.
	 * @param args
	 * @return the jobID
	 */
	public static int startJob(String[] args) {
		String JOBS_DIR = args[0];
		String jobType = args[2];
		
		String line = null;
		int jobID = 99;
		
		// get latest jobID
		BufferedReader in = null;
		try {
			in = new BufferedReader(new FileReader(JOBS_DIR + "/latestjob.id"));
			line = in.readLine();
		} catch (IOException e) {
			System.err.println("Error while reading the latestjob.id file.");
			e.printStackTrace();
		} finally {
			try {
				in.close();
			} catch (IOException e) {
				System.err.println("Error while trying to close the BufferedReader on letestjob.id file.");
				e.printStackTrace();
			}
		}
		jobID = Integer.parseInt(line);
		jobID++;
		
		// write latest jobID
		BufferedWriter out = null;
		line = String.valueOf(jobID);
		try {
			out = new BufferedWriter(new FileWriter(JOBS_DIR + "/latestjob.id"));
			out.write(line);
		} catch (IOException e) {
			System.err.println("Error while trying to write the latestjob.id file.");
			e.printStackTrace();
		} finally {
			try {
				out.close();
			} catch (IOException e) {
				System.err.println("Error while trying to close the BufferedWriter on latestjob.id file.");
				e.printStackTrace();
			}
		}
		
		// create JobFile $JOB_DIR/01_todo/$JOB_ID.job
		String todoFile = JOBS_DIR+"/01_todo/"+String.valueOf(jobID)+".job";
		try {
			out = new BufferedWriter(new FileWriter(todoFile));
			out.write("JOB_ID="+String.valueOf(jobID)+"\n");
			out.write("JOB_TYPE="+jobType+"\n");

			// TODO Write job specific arguments
			if (jobType.equalsIgnoreCase("tina")) {
				out.write("SEQUENCE="+args[3]);
			} else
			if (jobType.equalsIgnoreCase("gotoh")) {
				out.write("SEQUENCE_ONE="+args[3]+"\n");
				out.write("SEQUENCE_TWO="+args[4]+"\n");
				out.write("SHORT_MATRIX="+args[5]+"\n");
				out.write("MATRIX=\n"+args[6]);
			} else
			if (jobType.equalsIgnoreCase("coord")) {
				
			} else
			if (jobType.equalsIgnoreCase("123d")) {
				
			} else
			if (jobType.equalsIgnoreCase("kabsch")) {
				
			}
			
		} catch (IOException e) {
			System.err.println("Error while trying to write "+todoFile+".");
			e.printStackTrace();
		} finally {
			try {
				out.close();
			} catch (IOException e) {
				System.err.println("Error while trying to close the BufferedWriter on "+todoFile+".");
				e.printStackTrace();
			}
		}
		
		// start job (if not already working)
		File dir = new File(JOBS_DIR+"/02_working");
		if (dir.list().length < 1) {
			String[] arg2 = {JOBS_DIR};
			webservice.JobWorker.main(arg2);
		}
		
		return jobID;
		
	}
	
}
