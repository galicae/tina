/******************************************************************************
 * huberdp.oracles.ManualOracle.java                                          *
 * Contains the class ManualOracle which is an oracle for the RDP featuring   *
 * manual input.                                                              *
 *                                                                            *
 * This file is best read at line width 80 and tab width 4.                   *
 *                                                                   huberste *
 ******************************************************************************/
package huberdp.oracles;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.LinkedList;

import bioinfo.Sequence;
import bioinfo.alignment.SequenceAlignment;

import huberdp.Oracle;
import huberdp.PartialAlignment;
import huberdp.RDPProblem;

/**
 * ManualOracle is an Oracle for the HubeRDP which features manual input.
 * 
 * @author huberste
 * @lastchange 2013-02-14
 */
public class ManualOracle implements Oracle {

	@Override
	public LinkedList<PartialAlignment> findSimiliarSegments(
			RDPProblem problem, int m) {
		
		// allocate result
		LinkedList<PartialAlignment> results = new LinkedList<PartialAlignment>();

		// set Sequences for Oracle
		String[] rows = problem.getThreading().getRowsAsString();
		String template = "";
		String target = "";
		for (int i = problem.getProblemStart(); i < problem.getProblemEnd(); i++) {
			if (rows[0].charAt(i) != '-') {
				template += rows[0].charAt(i);
			}
			if (rows[1].charAt(i) != '-') {
				target += rows[1].charAt(i);
			}
		}
		
		Sequence templateSequence = new Sequence(problem.getThreading().getStructure().getID(), template);
		Sequence   targetSequence = new Sequence(problem.getThreading().getSequence().getId(), target);
		
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
			System.out.println("No further aligning is done.");
		} else {
			results.add(new PartialAlignment(problem, new SequenceAlignment(
					templateSequence, targetSequence, templateInput,
					targetInput, score)));
		}
		
		return null;
	}
}

/******************************************************************************
 * "A question that sometimes drives me hazy: * Am I or are the others crazy?" *
 * - Albert Einstein (1879 - 1955) *
 ******************************************************************************/
