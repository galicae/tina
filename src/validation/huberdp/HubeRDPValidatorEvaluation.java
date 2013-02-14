/******************************************************************************
 * validation.huberdp.HubeRDPValidatorEvaluation.java                         *
 *                                                                            *
 * This file is best read at line width 80 and tab width 4.                   *
 *                                                                   huberste *
 ******************************************************************************/
package validation.huberdp;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.Locale;

/**
 * HubeRDPValitadorEvaluation reads the output of an HubeRDPValidator and
 * generates an evaluation.
 * @author huberste
 * @lastchange 2013-02-14
 */
public class HubeRDPValidatorEvaluation {
	
	/**
	 * @param args
	 */
	public static void main(String[] args) {
		
		// For unified output
		Locale.setDefault(Locale.US);
		// DecimalFormat for formatting RMSD outputs
		DecimalFormat df = new DecimalFormat("0.0000");
		
		String inputfile = args[0];
		
		BufferedReader br = null;
		String line = null;
		String result = "";
		double huberdprmsd = 0.0;
		double gotohrmsd = 0.0;
		
		try {
			br = new BufferedReader(new FileReader(inputfile));
			while ((line = br.readLine()) != null) {
				if (line.startsWith(">>> ")) {
					if ((!line.startsWith(">>> HubeRDP:")) && (!line.startsWith(">>> Gotoh:"))) {
						String[] temp = line.split(" ");
						result = temp[1] + "\t" + temp[2];
					}
					
				} else if (line.startsWith("> HubeRDP RMSD:")) {
					huberdprmsd = Double.parseDouble(line.substring(16));
					result += "\t" + df.format(huberdprmsd);
				} else if (line.startsWith("> Gotoh RMSD:")) {
					gotohrmsd = Double.parseDouble(line.substring(14));
					result += "\t" + df.format(gotohrmsd);
					System.out.println(result);
					result = "";
				}
			}
		} catch (IOException e) {
			e.printStackTrace();
		} finally {
			try {
				br.close();
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}

	}

}

/******************************************************************************
 * "A question that sometimes drives me hazy:                                 *
 *  Am I or are the others crazy?"                                            *
 *     - Albert Einstein (1979 - 1966)                                        *
 ******************************************************************************/
