/******************************************************************************
 * bioinfo.energy.potential.hydrophobicity.HydrophobicityCalculator.java      *
 *                                                                            *
 * This class's main method uses the counts from                              *
 * bioinfo.energy.potential.hydrophobicity.DoBFreqCounter.main()              *
 * to calculate a hydrophobicity score based on the degree of burial (dob) of *
 * an AminoAcid and the frequencies of the dobs of all other AminoAcids.      *
 *                                                                            *
 * The idea is: w[i][j] = r[i][j]/(p[i]*q[j])                                 *
 * where r[i][j] is the relative frequency of amino acid [i] at dob[j],       *
 *       p[i] is the relative frequency of amino acid [i] and                 *
 *       q[j] is the relative frequency of dob[j]                             *
 *                                                                            *
 * @see Protein Threading by Recursive Dynamic Programming. JMB 290, 757-779  *
 *                                                                            *
 * This file is best read at line width 80 and tab width 4.                   *
 *                                                                   huberste *
 ******************************************************************************/
package bioinfo.energy.potential.hydrophobicity;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.Locale;

/**
 * @author huberste
 * @lastchange 2013-02-16
 */
public class HydrophobicityCalculator {

	public final static String usage = 
		"usage:\n" +
		"\tjava HydrobhobicityCalculator <infile> <outpath>\n\n" +
		"where <infile> is an output file of DoBFreqCounter, \n" +
		"<outpath> is a (writable) path to the folder that shall \n" +
		"contain the output files.";
	
	public final static double log2 = Math.log10(2);
	
	/**
	 * @param args
	 */
	public static void main(String[] args) {
		if (args.length < 2) {
			System.out.println(usage);
			System.exit(1);
		}
		
		// initialize arguments
		String infile = args[0];
		String outpath = args[1];
		int buckets = 0;
		
		// initialize important stuff
		long freq[][] = new long[26][]; 
		
		// read file
		BufferedReader br = null;
		try {
			br = new BufferedReader(new FileReader(infile));
			String line = null;
			int linenr = 0;
			while ((line = br.readLine()) != null) {
				String[] temp = line.split("\t");
				// buckets = buckets < temp.length-1 ? temp.length-1 : buckets;
				freq[linenr] = new long[temp.length-1];
				for (int i = 0; i < temp.length-1; i++) {
					freq[linenr][i] = Long.parseLong(temp[i+1]);
				}
				linenr++;
			}
		} catch (IOException e) {
			System.err.println("Error 62: problems reading the inputfile:");
			e.printStackTrace();
		} finally {
			try {
				br.close();
				br = null;
			} catch (IOException e) {
				System.err.println("Error 68: problems closing the inputfile:");
				e.printStackTrace();
			}
		}
		
		for(int size = 0; size < 1; size++) {
		
			buckets = 1024 / zweihoch(size);
			int sum = zweihoch(size);
			// crunch the numbers
			// sum up rows / cols
			long[] aacount = new long[26];
			long[] bucketcount = new long[buckets];
			
			long gescount = 0;
			for(int i = 0; i < 26; i++) {
				for (int j = 0; j < buckets; j++) {
					for(int k = 0; k < sum; k++) {
						aacount[i] += freq[i][(j*size)+k];
						bucketcount[j] += freq[i][(j*size)+k];
						gescount += freq[i][(j*size)+k];
					}
				}
			}
			
			// calculate r
			double[][] r = new double[26][buckets];
			for (int i = 0; i < 26; i++) {
				for (int j = 0; j < buckets; j++) {
					for(int k = 0; k < sum; k++) {
						r[i][j] = (double) freq[i][(j*size) +k] / (double) gescount;
					}
				}
			}
			
			// calculate p
			double[] p = new double[26];
			for (int i = 0; i < 26; i++) {
				p[i] = (double) aacount[i] / (double) gescount;
			}
			
			// calculate q
			double[] q = new double[buckets];
			for (int j = 0; j < buckets; j++) {
				q[j] = (double) bucketcount[j] / (double) gescount;
			}
		
		
		
			// TODO calculate w (s)
			double[][] w = new double[26][buckets];
			
			for(int i = 0; i < 26; i++) {
				for (int j = 0; j < buckets; j++) {
					double temp = r[i][j]/ (p[i]*q[j]);
					w[i][j] = Math.log10(temp) / log2;
				}
			}
			
			// TODO print out w
			print(w, outpath+"buckets.hmt");
		}
		
		
	}
	
	// TODO provide file structures and data structures for hydrophobicities
	
	/**
	 * 
	 * @param size
	 * @return
	 */
	private static int zweihoch(int size) {
		if (size == 0) {
			return 1;
		} else {
			return 2 * zweihoch(size-1);
		}
	}

	/**
	 * 
	 * @param w
	 */
	private static void print(double[][] w, String filepath) {
		
		// initialize important stuff	
		Locale.setDefault(Locale.US);
		DecimalFormat df = new DecimalFormat("0.000000");
		
		BufferedWriter bw = null;
		try {
			bw = new BufferedWriter(new FileWriter(filepath));
			bw.write("# Hydrophobicity matrix\n");
			bw.write("# degree of burial x amino acid\n");
			bw.write("# dob buckets: "+ w[0].length + "\n");
			for (int i = 0; i < w.length; i++) {
				char x = (char)(65+i);
				bw.write(x);
				for (int j = 0; j < w[i].length; j++) {
					bw.write("\t" + df.format(w[i][j]));
				}
				bw.write("\n");
			}
		} catch (IOException e) {
			System.err.println("Error 187: problems reading an outfile:");
			e.printStackTrace();
		} finally {
			try {
				bw.close();
			} catch (IOException e) {
				System.err.println("Error 193: problems closing an outfile:");
				e.printStackTrace();
			}
		}
	}

}

/******************************************************************************
 * "A question that sometimes drives me hazy:                                 *
 *  Am I or are the others crazy?"                                            *
 *     - Albert Einstein (1879 - 1955)                                        *
 ******************************************************************************/
