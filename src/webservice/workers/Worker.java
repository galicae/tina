package webservice.workers;

/**
 * @author huberste
 * @lastchange 2013-02-26
 */
public abstract class Worker {

	protected int JOB_ID;
	protected final String JOB_FILE;
	protected final String DONE_FILE;
	
	public Worker(String jobFile) {
		JOB_FILE = jobFile;
		JOB_ID=Integer.parseInt(JOB_FILE.substring(JOB_FILE.lastIndexOf('/')+1, JOB_FILE.lastIndexOf('.')));
		String temp = JOB_FILE.substring(0, JOB_FILE.lastIndexOf("/"));
		// TODO debugging: Check file names
		System.out.println("temp: "+temp);
		DONE_FILE = temp.substring(0,temp.lastIndexOf("/"))+"/03_done/"+String.valueOf(JOB_ID)+".id";
		System.out.println("DONE_FILE: "+DONE_FILE);
	}
	
	public abstract void work();
	
	protected abstract void readFile();
	
	protected abstract void writeResult();
	
}
