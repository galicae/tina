package bioinfo.energy.potential;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.OutputStreamWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Set;


import bioinfo.energy.potential.voroEval.VoroEvalDataPoint;
import bioinfo.energy.potential.voroEval.VoroEvalTree;
import bioinfo.energy.potential.voronoi.VoronoiData;
import bioinfo.proteins.AminoAcidName;
import bioinfo.proteins.DSSPEntry;
import bioinfo.proteins.DSSPFileReader;
import bioinfo.proteins.PDBEntry;

/**
 * based on voronoi surface neighborhood detection and dssp surface values
 * @author andreseitz
 * 
 */
public class DSSPSolventPotential extends AVoroPotential {

	/**
	 * potential contains the actual mean force potential [a][b][c] a contains
	 * aminoacid information (one letter code ascii - 65) of partner 1 b
	 * contains aminoacid information (one letter code ascii - 65) of partner 2
	 * c contains area of face between the two partners with the following
	 * classes smaller then 25,50,75,100,125,150,bigger then 150, where all
	 * values smaller then 1 have to be ignored
	 */
	private final int[] size = {7,26,26};
	private final String dsspFolder;
	private final double minContact;
	private final double mkT = -0.582d;
	private final double gridHullExtend;
	private final double gridDensity;
	private final double gridClash;

	public DSSPSolventPotential(String vorobin, String dsspFolder) {
		super(vorobin);
		initPotential(size);
		this.dsspFolder = dsspFolder;

		this.minContact = 1.0d;
		this.gridHullExtend = 2.0d;
		this.gridDensity = 3.0d;
		this.gridClash = 4.0d;
	}

	public DSSPSolventPotential(String vorobin, String tmpdir, String dsspFolder) {
		super(vorobin, tmpdir);
		initPotential(size);
		this.dsspFolder = dsspFolder;

		this.minContact = 1.0d;
		this.gridHullExtend = 2.0d;
		this.gridDensity = 3.0d;
		this.gridClash = 4.0d;

	}

	public DSSPSolventPotential(String vorobin, String dsspFolder, double minContact, double gridHullExtend, double gridDensity, double gridClash) {
		super(vorobin);
		initPotential(size);
		this.dsspFolder = dsspFolder;

		this.minContact = minContact;
		this.gridHullExtend = gridHullExtend;
		this.gridDensity = gridDensity;
		this.gridClash = gridClash;
	}

	public DSSPSolventPotential(String vorobin, String tmpdir, String dsspFolder, double minContact, double gridHullExtend, double gridDensity, double gridClash) {
		super(vorobin, tmpdir);
		initPotential(size);
		this.dsspFolder = dsspFolder;

		this.minContact = minContact;
		this.gridHullExtend = gridHullExtend;
		this.gridDensity = gridDensity;
		this.gridClash = gridClash;

	}

	// public DSSPSolventPotential(String filename,String vorobin, String
	// tmpdir){
	// super(vorobin, tmpdir);
	// this.potential = new double[26][26][7];
	// this.dsspFolder = null;
	//
	// this.minContact = 1.0d;
	// this.gridHullExtend = 2.0d;
	// this.gridDensity = 3.0d;
	// this.gridClash = 4.0d;
	//
	// this.readFromFile(filename);
	// }
	//
	// public DSSPSolventPotential(String filename,String vorobin){
	// super(vorobin);
	// this.potential = new double[26][26][7];
	// this.dsspFolder = null;
	//
	// this.minContact = 1.0d;
	// this.gridHullExtend = 2.0d;
	// this.gridDensity = 3.0d;
	// this.gridClash = 4.0d;
	//
	// this.readFromFile(filename);
	// }

	public void calculateEval(List<String> dsspIds, VoroEvalTree tree) {
		DSSPEntry dssp = null;
		DSSPFileReader reader = new DSSPFileReader(dsspFolder);
		VoronoiData data = null;
		Set<Integer> solventIds = null;
		Set<Integer> pepIds = null;
		int[] acc;
		HashMap<Integer, HashMap<Integer, Double>> faces;
		HashMap<Integer, Double> neighbors;
		boolean surfaceFlag = false;
		List<Integer> surfaceIds = new ArrayList<Integer>();
		double surfaceArea = 0.0d;
		AminoAcidName[] names;
		// int[] index;
		HashMap<Integer, Double> surfaces;

		for (String dsspId : dsspIds) {

			// read from files, can be changed to db later!
			dssp = reader.readFromFolderById(dsspId);
			names = dssp.getNames();
			// index = dssp.getResIndex();
			surfaces = new HashMap<Integer, Double>();

			data = this.prepareWithGrid(dssp, gridHullExtend, gridDensity, gridClash, minContact);
			acc = dssp.getAccesability();
			pepIds = data.getPepIds();
			solventIds = data.getOuterGridIds();
			faces = data.getFaces();

			// System.out.println("\tdecomp done ... "+faces.size()+" faces");

			for (int id1 : pepIds) {
				if (faces.get(id1) == null) {
					continue;
				}
				neighbors = faces.get(id1);
				surfaceFlag = false;
				surfaceArea = 0.0d;
				for (int id2 : neighbors.keySet()) {
					if (solventIds.contains(id2) && neighbors.get(id2) > minContact) {
						surfaceArea += neighbors.get(id2);
						surfaceFlag = true;
					}
				}
				if (surfaceFlag) {
					surfaces.put(id1, surfaceArea);
					surfaceIds.add(id1);
				}
			}

			for (int i = 0; i != pepIds.size(); i++) {
				if (surfaces.containsKey(i)) {
					tree.insertData(new VoroEvalDataPoint(names[i], dsspId, i, acc[i], surfaces.get(i)));
				} else {
					tree.insertData(new VoroEvalDataPoint(names[i], dsspId, i, acc[i], 0.0));
				}
			}
		}

	}

	public void calculateEval(List<String> dsspIds, String outLoc) {
		DSSPEntry dssp = null;
		DSSPFileReader reader = new DSSPFileReader(dsspFolder);
		VoronoiData data = null;
		Set<Integer> solventIds = null;
		Set<Integer> pepIds = null;
		int[] acc;
		HashMap<Integer, HashMap<Integer, Double>> faces;
		HashMap<Integer, Double> neighbors;
		boolean surfaceFlag = false;
		List<Integer> surfaceIds = new ArrayList<Integer>();
		double surfaceArea = 0.0d;
		AminoAcidName[] names;
		// int[] index;
		HashMap<Integer, Double> surfaces;
		BufferedWriter bw = null;

		for (String dsspId : dsspIds) {

			try {
				bw = new BufferedWriter(new OutputStreamWriter(new FileOutputStream(outLoc + dsspId + ".res")));
			} catch (FileNotFoundException e) {
				e.printStackTrace();
			}

			// read from files, can be changed to db later!
			dssp = reader.readFromFolderById(dsspId);
			names = dssp.getNames();
			// index = dssp.getResIndex();
			surfaces = new HashMap<Integer, Double>();

			data = this.prepareWithGrid(dssp, gridHullExtend, gridDensity, gridClash, minContact);
			acc = dssp.getAccesability();
			pepIds = data.getPepIds();
			solventIds = data.getOuterGridIds();
			faces = data.getFaces();

			// System.out.println("\tdecomp done ... "+faces.size()+" faces");

			for (int id1 : pepIds) {
				if (faces.get(id1) == null) {
					continue;
				}
				neighbors = faces.get(id1);
				surfaceFlag = false;
				surfaceArea = 0.0d;
				for (int id2 : neighbors.keySet()) {
					if (solventIds.contains(id2) && neighbors.get(id2) > minContact) {
						surfaceArea += neighbors.get(id2);
						surfaceFlag = true;
					}
				}
				if (surfaceFlag) {
					surfaces.put(id1, surfaceArea);
					surfaceIds.add(id1);
				}
			}
			try {
				for (int i = 0; i != pepIds.size(); i++) {
					if (surfaces.containsKey(i)) {
						bw.append(names[i] + "\t" + i + "\t" + acc[i] + "\t" + surfaces.get(i) + "\t" + Math.abs(acc[i] - surfaces.get(i)) + "\n");
					} else {
						bw.append(names[i] + "\t" + i + "\t" + acc[i] + "\t" + 0.0 + "\t" + Math.abs(acc[i]) + "\n");
					}
				}
				bw.flush();
				bw.close();
			} catch (Exception e) {
				e.printStackTrace();
			}
		}

	}

	public void calculateEval(List<String> dsspIds, VoroEvalTree tree, String outLoc) {
		DSSPEntry dssp = null;
		DSSPFileReader reader = new DSSPFileReader(dsspFolder);
		VoronoiData data = null;
		Set<Integer> solventIds = null;
		Set<Integer> pepIds = null;
		int[] acc;
		HashMap<Integer, HashMap<Integer, Double>> faces;
		HashMap<Integer, Double> neighbors;
		boolean surfaceFlag = false;
		List<Integer> surfaceIds = new ArrayList<Integer>();
		double surfaceArea = 0.0d;
		AminoAcidName[] names;
		// int[] index;
		HashMap<Integer, Double> surfaces;
		BufferedWriter bw = null;

		for (String dsspId : dsspIds) {

			try {
				bw = new BufferedWriter(new OutputStreamWriter(new FileOutputStream(outLoc + dsspId + ".res")));
			} catch (FileNotFoundException e) {
				e.printStackTrace();
			}

			// read from files, can be changed to db later!
			dssp = reader.readFromFolderById(dsspId);
			names = dssp.getNames();
			// index = dssp.getResIndex();
			surfaces = new HashMap<Integer, Double>();

			data = this.prepareWithGrid(dssp, gridHullExtend, gridDensity, gridClash, minContact);
			acc = dssp.getAccesability();
			pepIds = data.getPepIds();
			solventIds = data.getOuterGridIds();
			faces = data.getFaces();

			// System.out.println("\tdecomp done ... "+faces.size()+" faces");

			for (int id1 : pepIds) {
				if (faces.get(id1) == null) {
					continue;
				}
				neighbors = faces.get(id1);
				surfaceFlag = false;
				surfaceArea = 0.0d;
				for (int id2 : neighbors.keySet()) {
					if (solventIds.contains(id2) && neighbors.get(id2) > minContact) {
						surfaceArea += neighbors.get(id2);
						surfaceFlag = true;
					}
				}
				if (surfaceFlag) {
					surfaces.put(id1, surfaceArea);
					surfaceIds.add(id1);
				}
			}

			try {
				for (int i = 0; i != pepIds.size(); i++) {
					if (surfaces.containsKey(i)) {
						tree.insertData(new VoroEvalDataPoint(names[i], dsspId, i, acc[i], surfaces.get(i)));
						bw.append(names[i] + "\t" + i + "\t" + acc[i] + "\t" + surfaces.get(i) + "\t" + Math.abs(acc[i] - surfaces.get(i)) + "\t" + dsspId + "\n");

					} else {
						tree.insertData(new VoroEvalDataPoint(names[i], dsspId, i, acc[i], 0.0));
						bw.append(names[i] + "\t" + i + "\t" + acc[i] + "\t" + 0.0 + "\t" + Math.abs(acc[i]) + "\t" + dsspId + "\n");
					}
				}
				bw.flush();
				bw.close();
			} catch (Exception e) {
				e.printStackTrace();
			}
		}

	}

	@Override
	public void calculateFromDATA(List<String> dsspIds) {
		DSSPEntry dssp = null;
		DSSPFileReader reader = new DSSPFileReader(dsspFolder);
		VoronoiData data = null;
		Set<Integer> solventIds = null;
		Set<Integer> pepIds = null;
		int[] acc;
		HashMap<Integer, AminoAcidName> point;
		HashMap<Integer, HashMap<Integer, Double>> faces;
		HashMap<Integer, Double> neighbors;
		boolean surfaceFlag = false;
		List<Integer> surfaceIds = new ArrayList<Integer>();
		int tmp = 0;
		int p1;
		int p2;
		int[] stateCount = new int[7];
		int accessability = 0;

		for (String dsspId : dsspIds) {
			System.err.println(dsspId);
			// read from files, can be changed to db later!
			dssp = reader.readFromFolderById(dsspId);
			data = this.prepareWithGrid(dssp, gridHullExtend, gridDensity, gridClash, minContact);
			acc = dssp.getAccesability();
			pepIds = data.getPepIds();
			solventIds = data.getOuterGridIds();
			faces = data.getFaces();
			point = data.getAminos();

			for (int id1 : pepIds) {
				if (faces.get(id1) == null) {
					continue;
				}
				neighbors = faces.get(id1);
				surfaceFlag = false;
				for (int id2 : neighbors.keySet()) {
					if (solventIds.contains(id2) && neighbors.get(id2) > minContact) {
						surfaceFlag = true;
					}
				}
				if (surfaceFlag) {
					accessability = acc[id1];
					surfaceIds.add(id1);
					if (accessability > minContact) {
						for (int id2 : neighbors.keySet()) {
							if (!surfaceIds.contains(id2)) {
								continue;
							}
							tmp = 0;
							for (int i = 25; i <= 150; i += 25) {
								if (accessability <= i * 1.0d) {
									break;
								}
								tmp++;
							}

							p1 = point.get(id1).getOneLetterCode().charAt(0) - 65;
							p2 = point.get(id2).getOneLetterCode().charAt(0) - 65;
							int[] path1 = {tmp,p1,p2};
							potential.setValue(path1, potential.getByAddress(path1).getValue()+1.0d);
							int[] path2 = {tmp,p2,p1};
							potential.setValue(path2, potential.getByAddress(path2).getValue()+1.0d);

							stateCount[tmp]++;
						}
					}
				}
			}
		}

		int[] path = new int[3];
		for (int k = 0; k != 7; k++) {
			path[0] = k;
			for (int i = 0; i != 26; i++) {
				path[1] = i;
				for (int j = 0; j != 26; j++) {
					path[2] = j;
					potential.setValue(path, (mkT * Math.log((potential.getByAddress(path).getValue() + 1) / (stateCount[tmp]+1))));
				}
			}
		}

	}

	@Override
	public double scoreModel(PDBEntry model) {
		VoronoiData data = prepareWithGrid(model, gridHullExtend, gridDensity, gridClash, minContact);
		Set<Integer> solventIds = null;
		Set<Integer> pepIds = null;
		HashMap<Integer, AminoAcidName> amino;
		HashMap<Integer, HashMap<Integer, Double>> faces;
		HashMap<Integer, Double> neighbors;
		boolean surfaceFlag = false;
		List<Integer> surfaceIds = new ArrayList<Integer>();
		double surfaceArea = 0.0d;
		int tmp = 0;
		int p1;
		int p2;
		double score = 0.0d;

		voro.decomposite(data);
		faces = data.getFaces();
		amino = data.getAminos();
		pepIds = data.getPepIds();
		solventIds = data.getOuterGridIds();
		faces = data.getFaces();
		amino = data.getAminos();

		for (int id1 : pepIds) {
			if (faces.get(id1) == null) {
				continue;
			}
			neighbors = faces.get(id1);
			surfaceFlag = false;
			surfaceArea = 0.0d;
			for (int id2 : neighbors.keySet()) {
				if (solventIds.contains(id2) && neighbors.get(id2) > minContact) {
					surfaceArea += neighbors.get(id2);
					surfaceFlag = true;
				}
			}
			if (surfaceFlag) {
				surfaceIds.add(id1);
				if (surfaceArea > minContact) {
					for (int id2 : neighbors.keySet()) {
						if(!amino.containsKey(id2)){
							continue;
						}
						tmp = 0;
						for (int i = 25; i <= 150; i += 25) {
							if (surfaceArea <= i * 1.0d) {
								break;
							}
							tmp++;
						}
						p1 = amino.get(id1).getOneLetterCode().charAt(0) - 65;
						p2 = amino.get(id2).getOneLetterCode().charAt(0) - 65;
						int[] path = {tmp,p1,p2};
						score += potential.getByAddress(path).getValue();
					}
				}
			}
		}
		return score;
	}

	
	@Override
	public double[] getAminoScores(PDBEntry model){
		double[] scores = new double[model.length()];
		VoronoiData data = prepareWithGrid(model, gridHullExtend, gridDensity, gridClash, minContact);
		Set<Integer> solventIds = null;
		Set<Integer> pepIds = null;
		HashMap<Integer, AminoAcidName> amino;
		HashMap<Integer, HashMap<Integer, Double>> faces;
		HashMap<Integer, Double> neighbors;
		boolean surfaceFlag = false;
		List<Integer> surfaceIds = new ArrayList<Integer>();
		double surfaceArea = 0.0d;
		int tmp = 0;
		int p1;
		int p2;
		double score = 0.0d;

		voro.decomposite(data);
		faces = data.getFaces();
		amino = data.getAminos();
		pepIds = data.getPepIds();
		solventIds = data.getOuterGridIds();
		
		int pepIdsSize = pepIds.size();
		int l = 0;
		for (int id1 = 0; id1 != pepIdsSize; id1++) {
			if(!pepIds.contains(id1)){
				pepIdsSize++;
				continue;
			}
			if (faces.get(id1) == null) {
				continue;
			}
			neighbors = faces.get(id1);
			surfaceFlag = false;
			surfaceArea = 0.0d;
			for (int id2 : neighbors.keySet()) {
				if (solventIds.contains(id2) && neighbors.get(id2) > minContact) {
					surfaceArea += neighbors.get(id2);
					surfaceFlag = true;
				}
			}
			if (surfaceFlag) {
				surfaceIds.add(id1);
				if (surfaceArea > minContact) {
					score = 0.0d;
					for (int id2 : neighbors.keySet()) {
						if(!amino.containsKey(id2)){
							continue;
						}
						tmp = 0;
						for (int i = 25; i <= 150; i += 25) {
							if (surfaceArea <= i * 1.0d) {
								break;
							}
							tmp++;
						}
						p1 = amino.get(id1).getOneLetterCode().charAt(0) - 65;
						p2 = amino.get(id2).getOneLetterCode().charAt(0) - 65;
						int[] path = {tmp,p1,p2};
						score = potential.getByAddress(path).getValue();
					}
					scores[l] = score;
					l++;
				}
			}
		}
		return scores;
	}
	/**
	 * TEST main method
	 */
	public static void main(String[] args) {
		
		List<String> dsspIds = new ArrayList<String>();
		String dsspLoc = args[0];
		String vorobin = args[1];
		String outfile = args[2];
		String minContact = args[3];
		String gridExtend = args[4];
		String gridDensity = args[5];
		String gridClash = args[6];
		String outLocation = "/Users/andreseitz/Desktop/";

		File file = new File(dsspLoc);
		for (File f : file.listFiles()) {
			dsspIds.add(f.getName().substring(0, 7));
		}
		//dsspIds.add("1gaiA00");

		DSSPSolventPotential pot = new DSSPSolventPotential(vorobin, dsspLoc, Double.parseDouble(minContact), Double.parseDouble(gridExtend), Double.parseDouble(gridDensity), Double.parseDouble(gridClash));
		//VoroEvalTree eval = new VoroEvalTree(Double.parseDouble(minContact));
		//pot.calculateEval(dsspIds, eval, outLocation);
		//eval.printTree();
		pot.calculateFromDATA(dsspIds);
		pot.writeToFile(outfile);
		// System.out.println("done");

		// GridSolvensPotential pot = new
		// GridSolvensPotential(args[2],args[1],VoroPrepType.CC);
		// PDBFileReader reader = new PDBFileReader(pdbLoc);
		// PDBEntry pdb = reader.readFromFolderById("1j2xA00");
		// System.out.println(pot.scoreModel(pdb));
		//System.out.println((System.currentTimeMillis()-start)/1000+" sec");
	}

}
