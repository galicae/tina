/******************************************************************************
 * webservice.JobCreator                                                      *
 *                                                                            *
 * This file is best read at line width 80 and tab width 4.                   *
 ******************************************************************************/
package webservice;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;

/**
 * The JobCreator creates the diverse jobs that need to be calculated.
 * @author gobi_12_4
 * @date December 08, 2012
 * @version alpha
 */
public class JobCreator {
	
	/**
	 * This method creates a new jobfile in $JOBS_DIR.
	 * It also checks if there is already a worker running and calls a worker
	 * if not.
	 * @param args The first String must be the path to the $JOBS_DIR.
	 * The second String must be the type of the job ("gotoh", "123D", ...).
	 * All following Strings need to be the variables the job needs (e.g.
	 * sequences, costs, matrixtypes, ...)
	 * @return the jobID under which the job was created.
	 */
	public static int createJob(String[] args) {
		String JOBS_DIR=args[0];
		String jobType=args[1];
		
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
			out.write("JOB_TYPE="+jobType);

			// TODO Write job specific arguments
			
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
		
		File dir = new File(JOBS_DIR+"/02_working");
		// TODO foolproofing: check if dir is directory
		if (dir.list().length < 1) {
			// FIXME call worker. Do this as a new process!
			try {
				System.out.println("calling: java"+ " webservice.JobWorker " + JOBS_DIR);
//				ProcessBuilder pb = new ProcessBuilder("java", "webservice.JobWorker " + JOBS_DIR);
//				pb.start();
//				System.out.println(pb.command());
//				Runtime.getRuntime().exec("echo \"hallo\"");
				Runtime.getRuntime().exec("java webservice.JobWorker "+JOBS_DIR);
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}
		
		return jobID;
	}

}
