/******************************************************************************
 * ManualOracle is an oracle for the RDP featuring manual input.              *
 *                                                                            *
 * This file is best read at line width 80 and tab width 4.                   *
 *                                                                   huberste *
 ******************************************************************************/
package huberdp.oracles;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;

import bioinfo.Sequence;
import bioinfo.alignment.SequenceAlignment;

import huberdp.Oracle;

/**
 * @author huberste
 * @lastchange 2013-02-11
 *
 */
public class ManualOracle extends Oracle {

	/**
	 * this is the _real_ alignment step.
	 * @param templateSequence
	 * @param targetSequence
	 * @return the alignment of the two given sequences
	 */
	@Override
	protected SequenceAlignment align(Sequence templateSequence, Sequence targetSequence) {
		SequenceAlignment result = null;
		System.out.println("Please align the following Sequences:");
		System.out.println("Press ENTER without any text to break aligning here.");
		System.out.println("template: " + templateSequence.getSequenceAsString());
		System.out.println("target  : " + targetSequence.getSequenceAsString());
		String templateInput = null;
		String targetInput = null;
		int score = 0;
		try{
		    BufferedReader bufferRead = new BufferedReader(new InputStreamReader(System.in));
		    System.out.print("insert template: ");
		    templateInput = bufferRead.readLine();
		    if (! templateInput.equals("")) {
			    System.out.print("insert target:   ");
			    targetInput   = bufferRead.readLine();
			    System.out.print("insert score: ");
			    score = Integer.parseInt(bufferRead.readLine());
		    }
		}
		catch(IOException e)
		{
			e.printStackTrace();
		}
		
		if (templateInput.equals("")) {
//			result = null; // already happened in line 33
		} else {		
			result = new SequenceAlignment(
					templateSequence,
					targetSequence,
					templateInput,
					targetInput,
					score);
		}
			
		return result;
	}
	

}