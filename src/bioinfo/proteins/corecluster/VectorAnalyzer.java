package bioinfo.proteins.corecluster;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import bioinfo.proteins.PDBEntry;
import bioinfo.proteins.PDBFileReader;


public class VectorAnalyzer {

	private static final int X = 0, Y = 1, Z = 2;
	private PDBFileReader pdbFileReader = new PDBFileReader();

	public VectorAnalyzer(String pdbDir, String dsspFile, String outputDir)  {
		
		Map<String, char[]> dssp = CoreDefinition.parseDsspToThreeState(dsspFile);
		Map<String, VectorAnnotation> vectorMap = new HashMap<String, VectorAnnotation>();
		for (String s : dssp.keySet()) {
			vectorMap.put(s, getVectorAnnotationByDSSP(s, dssp.get(s), 4));
		}
		Pair<Double> helixStats = getStatistics(vectorMap, 'H');
		Pair<Double> sheetStats = getStatistics(vectorMap, 'E');
		// Pair<Double> helixStats = new Pair<Double>(2.199039948244969 ,
		// 0.17897073822129447);
		// Pair<Double> sheetStats = new Pair<Double>(1.1698974566256055 ,
		// 0.5671007013996289);

		System.out.println(helixStats);
		System.out.println(sheetStats);

		Map<String, VectorAnnotation> refinedVectorMap = new HashMap<String, VectorAnnotation>();
		for (String s : vectorMap.keySet())
			refinedVectorMap.put(s, refineVectorAnnotation(vectorMap.get(s), helixStats, sheetStats, 4));
	}

	private VectorAnnotation getFullDsspAnnotation(String id, char[] dssp) {
		VectorAnnotation result = new VectorAnnotation(id);
		List<CompPair<Integer>> dsspHelix = getSecondaryStructureElements(dssp, 'H', 0);
		for (CompPair<Integer> p : dsspHelix) {
			result.appendVector(new Curve(id, p, null, 'H'));
		}
		List<CompPair<Integer>> dsspSheet = getSecondaryStructureElements(dssp, 'E', 0);
		for (CompPair<Integer> p : dsspSheet) {
			result.appendVector(new Curve(id, p, null, 'E'));
		}
		return result;
	}

	private VectorAnnotation getVectorAnnotationByDSSP(String id, char[] dssp, int lengthThreshold)  {
		Set<PDBEntry> currentStructures = new HashSet<PDBEntry>();
		currentStructures = pdbFileReader.getPDBfromFileSplittedByChain(id);
		VectorAnnotation result = new VectorAnnotation(id);
		for (PDBEntry c : currentStructures) {
			result.appendVectors(getSecondaryStructureVectors(c, dssp, 'H', lengthThreshold));
			result.appendVectors(getSecondaryStructureVectors(c, dssp, 'E', lengthThreshold));
		}
		return result;
	}

	private VectorAnnotation refineVectorAnnotation(VectorAnnotation vectorAnnotation, Pair<Double> helixStats, Pair<Double> sheetStats, int lengthThreshold)  {
		Set<PDBEntry> currentStructures = new HashSet<PDBEntry>();
		currentStructures = pdbFileReader.getPDBfromFileSplittedByChain(vectorAnnotation.getId());
		VectorAnnotation result = new VectorAnnotation(vectorAnnotation.getId());
		for (PDBEntry c : currentStructures) {
			result.appendVectors(refineSecondaryStructureVectors(vectorAnnotation.getAllVectors(), c, helixStats.getX(), helixStats.getY(), lengthThreshold, 'H'));
			result.appendVectors(refineSecondaryStructureVectors(vectorAnnotation.getAllVectors(), c, sheetStats.getX(), sheetStats.getY(), lengthThreshold, 'E'));
		}
		return result;
	}

	private List<Curve> getSecondaryStructureVectors(PDBEntry c, char[] dssp, char sse, int lengthThreshold) {
		List<CompPair<Integer>> secondaryStructurePositions = getSecondaryStructureElements(dssp, sse, lengthThreshold);
		List<Curve> result = new ArrayList<Curve>();
		for (CompPair<Integer> p : secondaryStructurePositions) {
			result.add(VectorMath.getCurve(c, p, sse));
		}
		return result;
	}

	private List<Curve> refineSecondaryStructureVectors(List<Curve> curves, PDBEntry chain, double avg, double sd, int lengthThreshold, char sse) {
		List<Curve> result = new ArrayList<Curve>();
		for (Curve c : curves) {
			if (c.getType() != sse) {
				continue;
			}
			if (!VectorMath.isVectorNormal(c, chain, avg, sd)) {
				result.addAll(splitVector(c, chain, avg, sd, lengthThreshold, sse));
			} else {
				result.add(c);
			}
		}
		return result;
	}

	private List<Curve> splitVector(Curve curve, PDBEntry chain, double avg, double sd, int lengthThreshold, char sse) {
		// System.out.println(vector);
		List<InternalNode> internalNodes = new ArrayList<InternalNode>();
		for (int i = curve.getPosition().getX(); i < curve.getPosition().getY() - lengthThreshold; i++) {
			for (int j = i + lengthThreshold; j < curve.getPosition().getY(); j++) {
				Curve tempCurve = VectorMath.getCurve(chain, new CompPair<Integer>(i, j), sse);
				if (!VectorMath.isVectorNormal(tempCurve, chain, avg, sd)) {
					continue;
				}
				internalNodes.add(new InternalNode(tempCurve, VectorMath.getAverageDistance(tempCurve, chain)));
			}
		}
		List<Curve> result = new ArrayList<Curve>();
		internalNodes = getOptimalPath(internalNodes);
		for (InternalNode n : internalNodes) {
			result.add(n.getCurve());
		}
		if (result.isEmpty()) {
			result.add(curve);
		}
		return result;
	}

	public List<InternalNode> getOptimalPath(List<InternalNode> internalNodes) {
		internalNodes = buildGraph(internalNodes);
		// for (InternalNode n : nodes) {
		// System.out.println(n);
		// }
		// System.out.println();
		int bestLength = 0;

		InternalNode bestNode = null;
		for (InternalNode n : internalNodes) {
			if (n.getLength() > bestLength) {
				bestLength = n.getLength();
				bestNode = n;
			}
			if (n.getLength() == bestLength) {
				if (bestNode.getDistance() > n.getDistance()) {
					bestLength = n.getLength();
					bestNode = n;
				}
			}
		}
		InternalNode nextNode = bestNode;
		List<InternalNode> result = new ArrayList<InternalNode>();
		while (nextNode != null) {
			result.add(nextNode);
			nextNode = nextNode.getBestPrefix();
		}
		// for (InternalNode n : result) {
		// System.out.println(n);
		// }
		return result;
	}

	



	private Pair<Double> getStatistics(Map<String, VectorAnnotation> vectorMap, char sse)  {
		double mean = getMean(vectorMap, sse);
		double sd = getStandardDeviation(vectorMap, mean, sse);
		return new Pair<Double>(mean, sd);
	}

	private double getMean(Map<String, VectorAnnotation> vectorMap, char sse)  {
		double cumulation = 0.0;
		int counter = 0;
		for (String s : vectorMap.keySet()) {
			Set<PDBEntry> currentStructures = new HashSet<PDBEntry>();
			currentStructures = pdbFileReader.getPDBfromFileSplittedByChain(s);
			for (PDBEntry chain : currentStructures) {
				for (Curve curve : vectorMap.get(s).getAllVectors()) {
					if (curve.getType() == sse) {
						cumulation += VectorMath.getAverageDistance(curve, chain);
						counter++;
					}
				}
			}
		}
		return cumulation / counter;
	}

	private double getStandardDeviation(Map<String, VectorAnnotation> vectorMap, double avg, char sse)  {
		double cumulation = 0.0;
		int counter = 0;
		for (String s : vectorMap.keySet()) {
			Set<PDBEntry> currentStructures = new HashSet<PDBEntry>();
			currentStructures = pdbFileReader.getPDBfromFileSplittedByChain(s);
			for (PDBEntry chain : currentStructures) {
				for (Curve curve : vectorMap.get(s).getAllVectors()) {
					double distance = VectorMath.getAverageDistance(curve, chain);
					if (curve.getType() == sse) {
						cumulation += Math.pow(distance - avg, 2.0);
						counter++;
					}
				}
			}
		}
		return Math.sqrt(cumulation / (counter - 1));

	}

	private List<InternalNode> buildGraph(List<InternalNode> internalNodes) {
		for (InternalNode n : internalNodes) {
			for (InternalNode m : internalNodes) {
				if (m.getEnd() <= n.getStart()) {
					n.addPrefix(m);
				}
			}
		}
		return internalNodes;
	}

	private List<CompPair<Integer>> getSecondaryStructureElements(char[] dssp, char ss, int cutoff) {
		// identify all elements
		List<CompPair<Integer>> result = new ArrayList<CompPair<Integer>>();
		int start = 0, end = -1;
		for (int i = 0; i < dssp.length; i++) {
			if (dssp[i] == ss) {
				if (end == -1) {
					start = i;
				}
				end = i;
			}
			if (dssp[i] != ss) {
				if (end != -1) {
					result.add(new CompPair<Integer>(start, end));
					end = -1;
				}
			}
		}
		if (end > -1) {
			result.add(new CompPair<Integer>(start, end));
		}

		// Filtering all short elements lower or equal cutoff
		List<Integer> removeIndices = new ArrayList<Integer>();
		for (int i = 0; i < result.size(); i++) {
			if (result.get(i).getY() - result.get(i).getX() < cutoff) {
				removeIndices.add(i);
			}
		}
		Collections.reverse(removeIndices);
		for (int i : removeIndices) {
			result.remove(i);
		}
		return result;
	}

	public static void main(String args[]) {
		try {
			// new VectorAnalyzer("/Users/ike/Documents/test",
			// "/Users/ike/Documents/test/test.dssp",
			// "/Users/ike/Documents/test/out/");
			new VectorAnalyzer("data/PDB/", "data/core/all.dssp", "data/core/vectors/");
			// new VectorAnalyzer(new File(args[0]), args[1]);
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}
