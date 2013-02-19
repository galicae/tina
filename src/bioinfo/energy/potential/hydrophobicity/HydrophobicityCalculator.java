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
import java.io.FileReader;
import java.io.IOException;

/**
 * @author huberste
 * @lastchange 2013-02-18
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
	 * @throws IOException 
	 */
	public static void main(String[] args) throws IOException {
		if (args.length < 2) {
			System.out.println(usage);
			System.exit(1);
		}
		
		// initialize arguments
		String infile = args[0];
		String outpath = args[1];
		int aas = 26;
		int buckets = 1024;
		
		// initialize important stuff
		long freq[][] = new long[aas][buckets]; 
		
		// read file
		BufferedReader br = null;
		try {
			br = new BufferedReader(new FileReader(infile));
			String line = null;
			int linenr = 0;
			while ((line = br.readLine()) != null) {
				if (line.startsWith("#")) { // comment
					continue;
				}
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
				if(br != null) {
					br.close();
					br = null;
				}
			} catch (IOException e) {
				System.err.println("Error 68: problems closing the inputfile:");
				e.printStackTrace();
			}
		}
		
		// PSEUDOCOUNTS
		for (int aa = 0; aa < aas; aa++) {
			for (int bucket = 0; bucket < buckets; bucket ++) {
				if(aa != 1 && aa != 9 & aa != 14 && aa != 20 & aa != 23 & aa != 25) {
					freq[aa][bucket]++;
				}
			}
		}
		
		for(int size = 0; size < 7; size++) {
		
			int bucketsize = zweihoch(size);
			int numbuckets = buckets / bucketsize;
			// crunch the numbers
			// sum up rows / cols
			long[] aacount = new long[aas];
			long[] bucketcount = new long[numbuckets];
			
			long gescount = 0;
			for(int aa = 0; aa < aas; aa++) {
				for (int bucket = 0; bucket < numbuckets; bucket++) {
					for(int k = 0; k < bucketsize; k++) {
						aacount[aa] += freq[aa][(bucket*bucketsize)+k];
						bucketcount[bucket] += freq[aa][(bucket*bucketsize)+k];
						gescount += freq[aa][(bucket*bucketsize)+k];
					}
				}
			}
			
			// calculate r
			double[][] r = new double[aas][numbuckets];
			for (int aa = 0; aa < aas; aa++) {
				for (int bucket = 0; bucket < numbuckets; bucket++) {
					for(int k = 0; k < bucketsize; k++) {
						r[aa][bucket] = (double) freq[aa][(bucket*bucketsize) +k] / (double) gescount;
					}
				}
			}
			
			// calculate p
			double[] p = new double[aas];
			for (int aa = 0; aa < aas; aa++) {
				p[aa] = (double) aacount[aa] / (double) gescount;
			}
			
			// calculate q
			double[] q = new double[numbuckets];
			for (int bucket = 0; bucket < numbuckets; bucket++) {
				q[bucket] = (double) bucketcount[bucket] / (double) gescount;
			}
		
			// calculate w (s)
			double[][] w = new double[aas][numbuckets];
			
			for(int aa = 0; aa < aas; aa++) {
				for (int bucket = 0; bucket < numbuckets; bucket++) {
					double temp = r[aa][bucket]/ (p[aa]*q[bucket]);
					w[aa][bucket] = Math.log10(temp) / log2;
				}
			}
			
			// print out w
			HydrophobicityMatrix.writeFile(w, outpath+numbuckets+"buckets");
		}
	}
	
	/**
	 * 
	 * @param n the exponent
	 * @return 2^n
	 */
	private static int zweihoch(int n) {
		if (n == 0) {
			return 1;
		} else {
			return 2 * zweihoch(n-1);
		}
	}

}

/******************************************************************************
 * "A question that sometimes drives me hazy:                                 *
 *  Am I or are the others crazy?"                                            *
 *     - Albert Einstein (1879 - 1955)                                        *
 ******************************************************************************/
