package bioinfo.proteins.fr4gment;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.HashMap;
import java.util.LinkedList;

import bioinfo.Sequence;
import bioinfo.alignment.SequenceAlignment;
import bioinfo.proteins.PDBEntry;
import bioinfo.proteins.fragm3nt.ProteinFragment;
import bioinfo.superpos.PDBReduce;

/**
 * this class should group (almost) all functions needed to produce a complete
 * prediction given the cores of a protein
 * 
 * @author galicae
 * 
 */
public class LoopBaseline {
	private MultipleSuperposition ms;
	private SequenceAlignment input;
	private LinkedList<int[]> templCores;
	private LinkedList<int[]> queryCores;

	public LoopBaseline(SequenceAlignment input, String clusterDirectory,
			String guideFile) {
		this.input = input;
		System.out.println(input.toStringVerbose());
		this.ms = findMultipleSuperpos(clusterDirectory, guideFile);
	}

	/**
	 * this function searches for the actual cluster where the id of the
	 * template can be found, with the help pf a catalog where every protein in
	 * the cathscop set is mapped to its cathscop id (for fold recognition)
	 * 
	 * @param clusterDirectory
	 *            the file that contains all cluster folders
	 * @param guideFile
	 *            the file with the mapping
	 * @return a {@link MultipleSuperposition} from the correct cluster
	 */
	private MultipleSuperposition findMultipleSuperpos(String clusterDirectory,
			String guideFile) {
		if (!clusterDirectory.endsWith("/"))
			clusterDirectory += "/";
		// clusterDirectory: msp_cluster05
		// guideFile: cathscop.ids

		String pdbId = input.getComponent(0).getID();
		String fold = parseGuide(guideFile, pdbId);

		File cluster = new File(clusterDirectory + fold + "/");
		File[] clusters = cluster.listFiles();
		if (clusters.length == 0)
			return null;

		String clusterFile = "";
		for (File f : clusters) {
			if (grep(pdbId, f)) {
				clusterFile = f.getAbsolutePath();
			}
		}

		MultipleSuperposition ms = new MultipleSuperposition(clusterFile);
		return ms;
	}

	/**
	 * this function parses a file that contains all the ids from the cathscop
	 * inpairs and outpairs sets and their cath classification in order to find
	 * out to which fold the given (query) id belongs
	 * 
	 * @param guideFile
	 *            the directory of the ids file
	 * @param query
	 *            the id of which we want to determine the fold
	 * @return the fold of query
	 */
	private String parseGuide(String guideFile, String query) {
		try {
			BufferedReader re = new BufferedReader(new FileReader(guideFile));
			String cathId = "";
			String line = "";

			while ((line = re.readLine()) != null) {
				if (line.contains(query)) {
					cathId = line.split("\\s+")[1];
					break;
				}
			}
			re.close();

			String result = cathId.split("\\.")[1];
			return result;
		} catch (Exception e) {
			e.printStackTrace();
		}

		return null;
	}

	/**
	 * cheap wrapper for a grep call, since I have no real motivation to write
	 * another BufferedReader.
	 * 
	 * @param query
	 *            the word to look for
	 * @param f
	 *            the file to look into
	 * @return true if anything is found (dangerous, there might be errors...)
	 *         false else.
	 */
	private boolean grep(String query, File f) {
		try {
			BufferedReader re = new BufferedReader(new FileReader(f));
			String line = "";
			while ((line = re.readLine()) != null) {
				if (line.startsWith("MODEL")) {
					if (line.contains(query)) {
						re.close();
						return true;
					}
				}
			}
			re.close();
			return false;
		} catch (Exception e) {
			e.printStackTrace();
		}
		return false;
	}

	/**
	 * this function maps every loop between two core elements on the sequence
	 * of every protein in the superposition and saves the structure fragments
	 * from the superposition corresponding to the loop in a HashMap, the key
	 * being their starting point in the template sequence
	 * 
	 * @param ms
	 *            a multiple superposition where our template sequence/structure
	 *            is included
	 * @param input
	 *            the alignment of the core regions of the query structure and a
	 *            template. Template is at position 0.
	 * @return a HashMap as described above
	 */
	private HashMap<Integer, LinkedList<ProteinFragment>> findLoops(
			MultipleSuperposition ms, SequenceAlignment input) {
		ms.sort(input.getComponent(0).getID());
		Sequence seq1 = input.getComponent(0);

		double[][] coord = PDBReduce.reduceSinglePDB(ms.getStructures().get(0));
		ProteinFragment usedFrag = new ProteinFragment(seq1.getID(),
				seq1.getSequenceAsString(), coord, coord.length);
		templCores = calcCorePoints(input);
		deriveQueryCores();
		HashMap<Integer, LinkedList<ProteinFragment>> betweenCores = new HashMap<Integer, LinkedList<ProteinFragment>>();

		for (int j = 1; j < templCores.size(); j++) {
			int[] nextCore = templCores.get(j);
			int[] prevCore = templCores.get(j - 1);
			
			LinkedList<ProteinFragment> currentLoop = new LinkedList<ProteinFragment>();
			currentLoop.add(usedFrag.getPart(prevCore[1], nextCore[0] + 1));
			
			for (int i = 1; i < ms.getStructures().size(); i++) {
				PDBEntry pdbX = ms.getStructures().get(i);
				double[][] xCoord = PDBReduce.reduceSinglePDB(pdbX);
				ProteinFragment x = new ProteinFragment(pdbX.getID(), pdbX
						.getSequence().getSequenceAsString(), xCoord,
						pdbX.length());
				CoreSegmentGotoh got = new CoreSegmentGotoh(-1, -1, 0.5,
						usedFrag, x);

				Sequence xSequence = new Sequence(x.getID(), x.getSequence());

				got.align(seq1, xSequence);
				int[] xCore = got.traceback(prevCore[1], nextCore[0]);
				// take last point of preceding and first point of following
				// core as well - is important for later assembly
				xCore[0]--;
				xCore[1]++;
				currentLoop.add(x.getPart(xCore));
			}
			betweenCores.put(prevCore[1], currentLoop);
		}
		return betweenCores;
	}

	private void deriveQueryCores() {
		queryCores = new LinkedList<int[]>();
		int[][] aligned = input.getAlignedResidues();

		for (int[] core : templCores) {
			int[] qCore = new int[2];
			for (int i = 0; i < aligned[0].length; i++) {
				if (aligned[0][i] == core[0]) {
					qCore[0] = aligned[1][i];
					int diff = core[1] - core[0];
					qCore[1] = aligned[1][i] + diff;
					break;
				}
			}
			queryCores.add(qCore);
		}
	}

	/**
	 * this function scans the alignment and returns the start and end points
	 * (in the known sequence, sequence 0) of the aligned regions, thus
	 * identifying the core points
	 * 
	 * @param input
	 *            the alignment (seq0 is the one with known structure and
	 *            multiple superposition)
	 * @return the list of all starts and ends of core segments
	 */
	private LinkedList<int[]> calcCorePoints(SequenceAlignment input) {
		LinkedList<int[]> result = new LinkedList<int[]>();

		int[][] aligned = input.getAlignedResidues();
		int start = 0;
		int end = 0;

		for (int i = 1; i < aligned[0].length; i++) {
			int diff0 = aligned[0][i] - aligned[0][i - 1];
			int diff1 = aligned[1][i] - aligned[1][i - 1];

			if (diff0 == 1 && diff1 == 1) {

			} else {
				end = i - 1;
				int[] temp = { aligned[0][start], aligned[0][end] };
				result.add(temp);
				start = i;
			}
		}
		end = aligned[0].length - 1;
		int[] temp = { aligned[0][start], aligned[0][end] };
		result.add(temp);
		return result;
	}

	public MultipleSuperposition getMultipleSuperposition() {
		return ms;
	}

	public SequenceAlignment getTemplateAlignment() {
		return input;
	}

	public LinkedList<int[]> getCorePoints() {
		return templCores;
	}

	public ProteinFragment makePrediction() {
		HashMap<Integer, LinkedList<ProteinFragment>> loopFragments = findLoops(
				ms, input);

		// the coordinates. First fill the cores out.
		double[][] resultCoordinates = new double[input.getComponent(1)
				.getSequence().length][3];
		fillOutCores(resultCoordinates);

		// filter loop segments according to sequence length
		for (int i = 1; i < templCores.size(); i++) {
			int[] queryPrev = queryCores.get(i - 1);
			int[] queryNext = queryCores.get(i);

			int[] prevCore = templCores.get(i - 1);
			filterSequenceLength(loopFragments.get(prevCore[1]), queryNext[0]
					- queryPrev[1] - 1);
			if(loopFragments.get(prevCore[1]).size() == 0)
				loopFragments.remove(prevCore[1]);
			prevCore[0] = prevCore[1] = -1;
		}

		// filter loop segments according to euclidean distance cutoff
		for (int i = 1; i < templCores.size(); i++) {
			int[] prevCore = templCores.get(i - 1);
			int[] nextCore = templCores.get(i);
			filterDistanceCutoff(loopFragments.get(prevCore[1]), prevCore,
					nextCore, 2);
		}

		String seq = input.getComponent(1).getSequenceAsString();
		return new ProteinFragment("coreTest", seq, resultCoordinates,
				seq.length());
	}

	/**
	 * this function copies the coordinates of the aligned core segments in the
	 * appropriate places
	 * 
	 * @param coord
	 */
	private void fillOutCores(double[][] coord) {
		int[][] aligned = input.getAlignedResidues();
		ProteinFragment tempStructure = new ProteinFragment("template",
				PDBReduce.reduceSinglePDB(ms.getStructures().get(0)), 1);
		// the start of the core in sequence coordinates is already given for
		// the template; we need to find the aligned coordinate for the query
		for (int[] templatePos : templCores) {
			int tempStart = templatePos[0];
			int tempEnd = templatePos[1];

			// we don't need querEnd because there are no gaps; the length of
			// the copied section is determined by tempEnd
			int querStart = -1;

			for (int i = 0; i < aligned[0].length; i++) {
				if (aligned[0][i] == tempStart) {
					querStart = aligned[1][i];
					break;
				}
			}

			for (int i = 0; i <= tempEnd - tempStart; i++) {
				coord[querStart + i] = tempStructure.getAllResidues()[tempStart
						+ i];
			}
		}
	}

	/**
	 * this function goes through the loop fragments in a loop segment and
	 * normalizes according to sequence length. All possible fragments of right
	 * length are created and saved in the input list
	 * 
	 * @param loopSegment
	 *            the collection of loop fragments from the superposition that
	 *            corresponds to this locus
	 * @param loopLength
	 *            the desired sequence length
	 */
	private void filterSequenceLength(LinkedList<ProteinFragment> loopSegment,
			int loopLength) {
		LinkedList<ProteinFragment> rightLengthLoop = new LinkedList<ProteinFragment>();
		for (ProteinFragment f : loopSegment) {
			int curLoopLength = f.getAllResidues().length;
			if (curLoopLength >= loopLength) {
				for (int i = 0; i < curLoopLength - loopLength; i++) {
					rightLengthLoop.add(f.getPart(i, i + loopLength));
				}
			}
		}
		loopSegment.clear();
		for(ProteinFragment f: rightLengthLoop) {
			loopSegment.add(f);
		}
	}

	/**
	 * filters the distance cutoff of a set of fragments of the correct sequence
	 * length: they must (almost) match the space gap between two core elements
	 * 
	 * @param loopSegment
	 *            the collection of loop fragments from the superposition that
	 *            corresponds to this locus
	 * @param prevCore
	 *            the last core element before the loop
	 * @param nextCore
	 *            the core element following the loop
	 * @param cutoff
	 *            the distance cutoff
	 */
	private void filterDistanceCutoff(LinkedList<ProteinFragment> loopSegment,
			int[] prevCore, int[] nextCore, double cutoff) {
		if(prevCore[0] < 0)
			return;
		// calculate the distance the core coordinates imply
		double[][] coord = PDBReduce.reduceSinglePDB(ms.getStructures().get(0));
		double[] start = coord[prevCore[1]];
		double[] end = coord[nextCore[0]];

		double dist = euclideanDistance(start, end);
		int fragLength = loopSegment.get(0).getAllResidues().length;

		LinkedList<ProteinFragment> rightDist = new LinkedList<ProteinFragment>();
		for (ProteinFragment f : loopSegment) {
			double curDist = euclideanDistance(f.getResidue(0),
					f.getResidue(fragLength - 1));
			if (isInCutoff(dist, curDist, cutoff))
				rightDist.add(f);
		}
	}

	/**
	 * method for euclidean distance calculation between two vectors of same
	 * dimensions
	 * 
	 * @param x
	 * @param y
	 * @return
	 */
	private double euclideanDistance(double[] x, double[] y) {
		double result = 0;
		for (int i = 0; i < x.length; i++) {
			result += Math.pow(x[i] - y[i], 2);
		}
		return Math.sqrt(result);
	}

	/**
	 * returns true if the difference of a and b is below cutoff
	 * 
	 * @param a
	 * @param b
	 * @param cutoff
	 * @return
	 */
	private boolean isInCutoff(double a, double b, double cutoff) {
		if (Math.abs(a - b) < cutoff)
			return true;
		else
			return false;
	}
}
