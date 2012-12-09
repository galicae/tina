/******************************************************************************
 * webservice.workers.TinaWorker                                              *
 *                                                                            *
 * This file is best read at line width 80 and tab width 4.                   *
 ******************************************************************************/
package webservice.workers;

/**
 * TinaWorker is the worker that works on an TINA job.
 * @author huberste
 * @date December 09, 2012
 * @version alpha
 */
public class TinaWorker extends Worker {
	
	/**
	 * 
	 * @param jobfile
	 */
	public TinaWorker(String jobfile) {
		super(jobfile);
	}
	
	/**
	 * 
	 */
	public void work() {
		// TODO save results in $JOB_DIR/03_done/$JOB_ID.job
		
		// TODO rm $JOB_DIR/02_working/$JOB_ID.job
	}

	@Override
	protected void readFile() {
		// TODO Auto-generated method stub
		
	}

}
