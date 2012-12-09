/**
 * 
 */
package webservice.workers;

/**
 * @author huberste
 *
 */
public abstract class Worker {

	protected int JOB_ID;
	protected final String JOB_FILE;
	
	public Worker(String jobFile) {
		JOB_FILE = jobFile;
		JOB_ID=Integer.parseInt(JOB_FILE.substring(JOB_FILE.lastIndexOf('/'), JOB_FILE.lastIndexOf('.')));
	}
	
	public abstract void work();
	
	protected abstract void readFile();
	
}
