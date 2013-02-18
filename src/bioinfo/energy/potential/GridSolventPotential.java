package bioinfo.energy.potential;


import java.io.File;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Set;
import bioinfo.energy.potential.voronoi.VoronoiData;
import bioinfo.proteins.AminoAcidName;
import bioinfo.proteins.PDBEntry;
import bioinfo.proteins.PDBFileReader;

public class GridSolventPotential extends AVoroPotential {

	/**
	 * potential contains the actual mean force potential [a][b][c] a contains
	 * aminoacid information (one letter code ascii - 65) of partner 1 b
	 * contains aminoacid information (one letter code ascii - 65) of partner 2
	 * c contains area of face between the two partners with the following
	 * classes smaller then 25,50,75,100,125,150,bigger then 150, where all
	 * values smaller then 1 have to be ignored
	 */
	private final String pdbFolder;
	private final int[] size = {7,26,26};
	
	private final double minContact;
	private final double gridHullExtend;
	private final double gridDensity;
	private final double gridClash;
	private final double mkT = -0.582d;

	public GridSolventPotential(String vorobin, String pdbFolder, double minContact, double gridHullExtend, double gridDensity, double gridClash) {
		super(vorobin);
		initPotential(size);
		this.pdbFolder = pdbFolder;

		this.minContact = minContact;
		this.gridHullExtend = gridHullExtend;
		this.gridDensity = gridDensity;
		this.gridClash = gridClash;
	}

	public GridSolventPotential(String vorobin, String tmpdir, String pdbFolder, double minContact, double gridHullExtend, double gridDensity, double gridClash) {
		super(vorobin, tmpdir);
		initPotential(size);
		this.pdbFolder = pdbFolder;

		this.minContact = minContact;
		this.gridHullExtend = gridHullExtend;
		this.gridDensity = gridDensity;
		this.gridClash = gridClash;

	}

	public GridSolventPotential(String filename, String vorobin, String tmpdir) {
		super(vorobin, tmpdir);
		initPotential(size);
		this.pdbFolder = null;

		this.minContact = 1.0d;
		this.gridHullExtend = 1.0d;
		this.gridDensity = 1.0d;
		this.gridClash = 1.0d;
	}

	public GridSolventPotential(String filename, String vorobin) {
		super(vorobin);
		initPotential(size);
		this.pdbFolder = null;

		minContact = 1.0d;
		gridHullExtend = 1.0d;
		gridDensity = 1.0d;
		gridClash = 1.0d;
	}

	public void calculateFromDATA(List<String> pdbIds) {
		PDBEntry pdb = null;
		PDBFileReader reader = new PDBFileReader(pdbFolder);
		VoronoiData data = null;
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
		int[] stateCount = new int[7];

		for (String pdbId : pdbIds) {

			// read from files, can be changed to db later!
			pdb = reader.readFromFolderById(pdbId);
			data = this.prepareWithGrid(pdb, gridHullExtend, gridDensity, gridClash, minContact);

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
					potential.setValue(path, (mkT * Math.log(((potential.getByAddress(path).getValue() + 1)) / stateCount[k])));
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

		List<String> pdbIds = new ArrayList<String>();
		String pdbLoc = args[0];
		String vorobin = args[1];
		String outfile = args[2];
		String minContact = args[3];
		String gridExtend = args[4];
		String gridDensity = args[5];
		String gridClash = args[6];

		File file = new File(pdbLoc);
		for (File f : file.listFiles()) {
			pdbIds.add(f.getName().substring(0, 7));
		}
		//pdbIds.add("1j2xA00");


		GridSolventPotential pot = new GridSolventPotential(vorobin, pdbLoc, Double.parseDouble(minContact), Double.parseDouble(gridExtend), Double.parseDouble(gridDensity), Double.parseDouble(gridClash));
		pot.calculateFromDATA(pdbIds);
		pot.writeToFile(outfile);
		//System.out.println("done");

		// GridSolvensPotential pot = new
		// GridSolvensPotential(args[2],args[1],VoroPrepType.CC);
		// PDBFileReader reader = new PDBFileReader(pdbLoc);
		// PDBEntry pdb = reader.readFromFolderById("1j2xA00");
		// System.out.println(pot.scoreModel(pdb));

	}

}
