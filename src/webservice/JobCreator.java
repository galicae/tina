/******************************************************************************
 * webservice.JobCreator                                                      *
 *                                                                            *
 * This file is best read at line width 80 and tab width 4.                   *
 ******************************************************************************/
package webservice;

import java.io.BufferedInputStream;
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
	 * Second String must be the $WEBAPP_HOME.
	 * Third String must be the type of the job ("gotoh", "123D", ...).
	 * All following Strings need to be the variables the job needs (e.g.
	 * sequences, costs, matrixtypes, ...)
	 * @return the jobID under which the job was created.
	 */
	public static int createJob(String[] args) {
		String JOBS_DIR = args[0];
		String WEBAPP_HOME = args[1];
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
			try {
				String workingDir = WEBAPP_HOME + "/WEB-INF/classes";
				String arg1 = "webservice.JobWorker";
				String arg2 = JOBS_DIR;

				ProcessBuilder pb = new ProcessBuilder("java", arg1, arg2);
				pb.directory(new File(workingDir));
				
				// DONE debugging: what command was called?
//					System.out.println("New process will be started now: "+pb.command());
				// end debugging
				
				pb.start();
//				Process proc = pb.start();
//
//				// TODO debugging: Give error stream! GIVE!
//					BufferedInputStream err = new BufferedInputStream( proc.getErrorStream() );
//					BufferedInputStream outstr = new BufferedInputStream( proc.getInputStream() );
//					// print output
//					byte[] buf = new byte[1024];
//					int nr = outstr.read(buf);
//					while (nr != -1)
//					{
//						System.out.write(buf, 0, nr);
//						nr = outstr.read(buf);
//					}
//					nr = err.read(buf);
//					while (nr != -1)
//					{
//						System.err.write(buf, 0, nr);
//						nr = err.read(buf);
//					}
				// end debugging
			} catch (Exception e) {
				e.printStackTrace();
			}
		}
		
		return jobID;
	}

}
