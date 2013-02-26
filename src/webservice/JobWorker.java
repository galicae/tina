/******************************************************************************
 * webservice.JobWorker                                                       *
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

import webservice.workers.*;

/**
 * The JobWorker works on created jobs.
 * @author gobi_12_4
 * @date December 08, 2012
 * @version alpha
 */
public class JobWorker {

	// TODO paralellize workers!
	// TODO make workers work in a grid!
	
	private static final int WORKER_LIMIT = 1;
	
	/** 
	 * function to call a worker.
	 * Workers should not be called directly but only through this main method.
	 * @param args the first argument must be $JOB_DIR. No arguments must follow.
	 */
	public static void main (String[] args) {
//		System.out.println("Line 36");
		
		String JOBS_DIR = args[0];
		JobWorker[] workers = new JobWorker[WORKER_LIMIT];
		int workerCount = 0;
		int jobID = -1;
		
		File dir = new File(JOBS_DIR+"/01_todo");
		// while worker limit not reached &&   there is work to do
		while ((workerCount < WORKER_LIMIT) && (dir.list().length > 0)) {
			String jobFile = dir.list()[0];
			jobID = Integer.parseInt(jobFile.split("\\.")[0]);
			
			workers[workerCount] = new JobWorker(JOBS_DIR, jobID);
			workers[workerCount].work();
		}

		
	}
	
	private final String JOB_DIR;
	private final int JOB_ID;
	
	/**
	 * creates a new worker
	 * @param jobdir
	 * @param jobid
	 */
	public JobWorker(String jobdir, int jobid) {
		JOB_DIR = jobdir;
		JOB_ID = jobid;
	}
	
	/**
	 * starts the worker
	 * here all the magic happens.
	 */
	private void work() {
		
		String todoFile = JOB_DIR + "/01_todo/"+JOB_ID+".job";
		String workingFile = JOB_DIR + "/02_working/"+JOB_ID+".job";
//		String doneFile = JOB_DIR + "/03_done/"+JOB_ID+".job";
		
		// mv $JOB_DIR/01_todo/$JOB_ID.job $JOB_DIR/02_working/$JOB_ID.job
		// && read $JOB_TYPE
		String jobType = null; 
		BufferedReader from = null;
		BufferedWriter to = null;
		
		// copy file to 02_working
		String line = null; 
		try {
			from = new BufferedReader(new FileReader(todoFile));
			to = new BufferedWriter(new FileWriter(workingFile));
			while((line = from.readLine()) != null) {
				to.write(line+"\n");
				// read $JOB_TYPE
				if (line.startsWith("JOB_TYPE=")) {
					jobType=line.substring(9);
				}
				
			}
		} catch (IOException e) {
			System.err.println("Error while trying to copy "+todoFile+" to "+workingFile+".");
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
		// delete jobfile in 01_todo
		new File(todoFile).delete();
		
		Worker worker = null;
		// TODO add other workers
		if (jobType.equalsIgnoreCase("tina")) {
			worker = new TinaWorker(workingFile);
		} else if (jobType.equalsIgnoreCase("gotoh")) {
			worker = new GotohWorker(workingFile);
		} else if (jobType.equalsIgnoreCase("coord")) {
			worker = new CoordWorker(workingFile);
		} else if (jobType.equalsIgnoreCase("123d")) {
			worker = new D123Worker(workingFile);
		} else if (jobType.equalsIgnoreCase("kabsch")) {
			worker = new KabschWorker(workingFile);
		} else if (jobType.equalsIgnoreCase("huberdp")) {
			worker = new HubeRDPWorker(workingFile);
		} else if (jobType.equalsIgnoreCase("kerbsch")) {
			worker = new KerbschWorker(workingFile);
		} else if (jobType.equalsIgnoreCase("muster")) {
			worker = new MusterWorker(workingFile);
		} else if (jobType.equalsIgnoreCase("fragm3nt")) {
			worker = new Fragm3ntWorker(workingFile);
		} else if (jobType.equalsIgnoreCase("baseline")) {
			worker = new BaselineWorker(workingFile);
		}
		
		// Do the actual work
		worker.work();
		
	}
	
	/**
	 * Copies a File in a very inefficient way
	 * @param fromFile
	 * @param toFile
	 * @return false if failed, else true
	 */
	private static boolean copy(String fromFile, String toFile) {
		BufferedReader from = null;
		BufferedWriter to = null;
		
		String line = null; 
		try {
			from = new BufferedReader(new FileReader(fromFile));
			to = new BufferedWriter(new FileWriter(toFile));
			while((line = from.readLine()) != null) {
				to.write(line+"\n");
			}
		} catch (IOException e) {
			System.err.println("Error while trying to copy "+fromFile+" to "+toFile+".");
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
		
		// TODO return false is something went wrong!
		return true;
	}
	
	/**
	 * Moves a file in a very, very inefficient way
	 * @param fromFile
	 * @param toFile
	 * @return false if failed, else true
	 */
	static boolean move(String fromFile, String toFile) {
		boolean result = (copy(fromFile, toFile) && (new File(fromFile).delete()));
		return result;
	}
	
}
