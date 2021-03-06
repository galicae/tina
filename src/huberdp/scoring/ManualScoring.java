/******************************************************************************
 * huberdp.scoring.ManualScoring.java                                         *
 *                                                                            *
 * Contains the class ManualScoring that provides a function to manually      *
 * score OR nodes.                                                            *
 *                                                                            *
 * This file is best read at line width 80 and tab width 4.                   *
 *                                                                   huberste *
 ******************************************************************************/
package huberdp.scoring;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;

import bioinfo.alignment.Threading;
import huberdp.Scoring;

/**
 * ManualScoring provides a function to manually score OR nodes.
 * @author huberste
 * @lastchange 2013-02-14
 */
public class ManualScoring implements Scoring {

	/* (non-Javadoc)
	 * @see huberdp.Scoring#score(huberdp.RDPSolutionTreeOrNode)
	 */
	@Override
	public double score(Threading threading) {
		System.out.println("Please score the following Threading.");
		System.out.println("Press ENTER to set default score (0.0).");
		System.out.println(threading);
		double score = 0.0;
		try{
		    BufferedReader bufferRead = new BufferedReader(new InputStreamReader(System.in));
		    System.out.print("insert score: ");
		    String temp = bufferRead.readLine();
		    if (!temp.equals("")) {
			    score = Double.parseDouble(temp);
		    }
		}
		catch(IOException e)
		{
			e.printStackTrace();
		}
		return score;
	}

	@Override
	public double getScore(Threading t, int m, int n) {
		// TODO Auto-generated method stub
		return 0;
	}

	@Override
	public double getInsertionScore(Threading t, int n) {
		// TODO Auto-generated method stub
		return 0;
	}

	@Override
	public double getDeletionScore(Threading t, int m) {
		// TODO Auto-generated method stub
		return 0;
	}

}

/******************************************************************************
 * "A question that sometimes drives me hazy:                                 *
 *  Am I or are the others crazy?"                                            *
 *     - Albert Einstein (1879 - 1955)                                        *
 ******************************************************************************/
