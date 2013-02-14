package bioinfo.energy.potential;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.InputStreamReader;
import java.io.OutputStreamWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import bioinfo.energy.potential.voronoi.VoroPPWrap;
import bioinfo.energy.potential.voronoi.VoroPrepType;
import bioinfo.energy.potential.voronoi.VoronoiData;
import bioinfo.proteins.AminoAcidName;
import bioinfo.proteins.PDBEntry;
import bioinfo.proteins.PDBFileReader;

public class GridSolventPotential extends AVoroPotential{

	/**
	 * potential contains the actual mean force potential
	 * [a][b][c]
	 * a contains aminoacid information (one letter code ascii - 65) of partner 1
	 * b contains aminoacid information (one letter code ascii - 65) of partner 2
	 * c contains area of face between the two partners with the following classes
	 * smaller then 25,50,75,100,125,150,bigger then 150, where all values smaller then 1 have to be ignored
	 */
	private final String pdbFolder;
	private final VoroPrepType type;
	
	private final double minContact;
	private final double gridHullExtend;
	private final double gridDensity;
	private final double gridClash;
	private final double mkT = -0.582d;
	
	private final String[] mappingKeys = {"aminoacid1","aminoacid2","faceArea"};
	
	public GridSolventPotential(String vorobin, String pdbFolder, List<String> pdbIds, VoroPrepType type){
		super(vorobin);
		this.potential = new double[26][26][7];
		this.pdbFolder = pdbFolder;
		this.type = type;
		
		minContact = 1.0d;
		gridHullExtend = 2.0d;
		gridDensity = 3.0d;
		gridClash = 4.0d;
		
		calculateFromDATA(pdbIds);
		
	}
	
	public GridSolventPotential(String vorobin, String tmpdir, String pdbFolder, List<String> pdbIds, VoroPrepType type){
		super(vorobin,tmpdir);
		this.potential = new double[26][26][7];
		this.pdbFolder = pdbFolder;
		this.type = type;
		
		minContact = 1.0d;
		gridHullExtend = 2.0d;
		gridDensity = 3.0d;
		gridClash = 4.0d;
		
		calculateFromDATA(pdbIds);

	}
	
	public GridSolventPotential(String vorobin, String pdbFolder, List<String> pdbIds, VoroPrepType type, double minContact, double gridHullExtend, double gridDensity, double gridClash){
		super(vorobin);
		this.potential = new double[26][26][7];
		this.pdbFolder = pdbFolder;
		this.type = type;
		
		this.minContact = minContact;
		this.gridHullExtend = gridHullExtend;
		this.gridDensity = gridDensity;
		this.gridClash = gridClash;
		
		calculateFromDATA(pdbIds);
		
	}
	
	public GridSolventPotential(String vorobin, String tmpdir, String pdbFolder, List<String> pdbIds, VoroPrepType type, double minContact, double gridHullExtend, double gridDensity, double gridClash){
		super(vorobin, tmpdir);
		this.potential = new double[26][26][7];
		this.pdbFolder = pdbFolder;
		this.type = type;

		this.minContact = minContact;
		this.gridHullExtend = gridHullExtend;
		this.gridDensity = gridDensity;
		this.gridClash = gridClash;
		
		calculateFromDATA(pdbIds);

	}
	
	public GridSolventPotential(String filename,String vorobin, String tmpdir,VoroPrepType type){
		super(vorobin,tmpdir);
		this.potential = new double[26][26][7];
		this.type=type;
		this.pdbFolder=null;
		
		this.minContact = 1.0d;
		this.gridHullExtend = 2.0d;
		this.gridDensity = 3.0d;
		this.gridClash = 4.0d;
		
		readFromFile(filename);
	}
	
	public GridSolventPotential(String filename,String vorobin,VoroPrepType type){
		super(vorobin);
		this.potential = new double[26][26][7];
		this.type=type;
		this.pdbFolder=null;
		
		minContact = 1.0d;
		gridHullExtend = 2.0d;
		gridDensity = 3.0d;
		gridClash = 4.0d;
		
		readFromFile(filename);
	}

	public void calculateFromDATA(List<String> pdbIds) {
		PDBEntry pdb = null;
		PDBFileReader reader = new PDBFileReader(pdbFolder);
		VoronoiData data = null;
		Set<Integer> solventIds = null;
		Set<Integer> pepIds = null;
		HashMap<Integer, AminoAcidName> amino;
		HashMap<Integer, HashMap<Integer,Double>> faces;
		HashMap<Integer,Double> neighbors;
		boolean surfaceFlag = false;
		List<Integer> surfaceIds = new ArrayList<Integer>();
		int tmp = 0;
		int p1;
		int p2;
		int[] aminoCount = new int[26];
		double surfaceArea;

		
		
		for(String pdbId: pdbIds){
			
			System.out.println(">"+pdbId);
			
			//read from files, can be changed to db later!
			pdb = reader.readFromFolderById(pdbId);
			data = this.prepareWithGrid(pdb, type, gridHullExtend, gridDensity, gridClash, minContact);
			
			pepIds = data.getPepIds();
			solventIds = data.getOuterGridIds();
			faces = data.getFaces();
			amino = data.getAminos();
			
			System.out.println("\tdecomp done ... "+faces.size()+" faces");
			
			for(int id1: pepIds){
				if(faces.get(id1) == null){
					continue;
				}
				neighbors = faces.get(id1);
				surfaceFlag = false;
				surfaceArea = 0.0d;
				for(int id2 : neighbors.keySet()){
					if(solventIds.contains(id2) && neighbors.get(id2) > minContact){
						surfaceArea += neighbors.get(id2);
						surfaceFlag = true;
					}
				}
				if(surfaceFlag){
					//System.err.println(surfaceArea);
					surfaceIds.add(id1);
					if(surfaceArea > minContact){
						for(int id2 : neighbors.keySet()){
							tmp = 0;
							for(int i = 25; i <= 150; i += 25){
								if(surfaceArea <= i*1.0d){
									break;
								}
								tmp++;
							}
								p1 = amino.get(id1).getOneLetterCode().charAt(0)-65;
								p2 = amino.get(id2).getOneLetterCode().charAt(0)-65;
								potential[p1][p2][tmp]++;
								potential[p2][p1][tmp]++;
								aminoCount[p1]++;
							}
						}
					}
				}
			}
			
		
		
		for(int i = 0; i != 26; i++){
			for(int j = 0; j != 26; j++){
				for(int k = 0; k != 7; k++){
					if(i == j){
						potential[i][j][k] = mkT*Math.log(((potential[i][j][k]+1)/2)/(((aminoCount[i]+aminoCount[j])/2)+1));
					}else{
						potential[i][j][k] = mkT*Math.log((potential[i][j][k]+1)/(((aminoCount[i]+aminoCount[j])/2)+1));
					}
				}
			}
		}
		
	}

	/**
	 * function to access single energy values from potential
	 * dont use as batch accessing method - might be very slow due to parameter passing as map
	 * key values for map can be looked up in potentials class commi
	 * 21.01.1013 required keys are aminoacid1,aminoacid2 and faceArea
	 */
	@Override
	public double getEnergyValue(Map<Object, Object> mapping) {
		String missing_keys = "";
		for(String key : this.mappingKeys){
			if(!mapping.containsKey(key)){
				missing_keys += key+"\t";
			}
		}
		if(missing_keys.length() > 0){
			System.err.println("key(s) missing: "+missing_keys+" in energy value request!");
			return 0.0d;
		}
		
		char aa1 = (Character)mapping.get("aminoacid1");
		char aa2 = (Character)mapping.get("aminoacid2");
		double area = (Double)mapping.get("faceArea");
		int k = 0;
		for(int i = 25; i <= 150; i += 25){
			if(area <= i*1.0d){
				break;
			}
			k++;
		}
		
		int i = aa1-65;
		int j = aa2-65;
		
		return potential[i][j][k];
		
	}

	@Override
	public String[] getMapKeys() {
		return this.mappingKeys;
	}

	@Override
	public double scoreModel(PDBEntry model) {
		VoronoiData data = prepareWithGrid(model, type, gridHullExtend, gridDensity, gridClash, minContact);
		Set<Integer> pepIds = data.getPepIds();
		Set<Integer> gridIds = data.getOuterGridIds();
		HashMap<Integer, AminoAcidName> amino;
		HashMap<Integer, HashMap<Integer,Double>> faces;
		HashMap<Integer,Double> neighbors;
		boolean surfaceFlag = false;
		List<Integer> surfaceIds = new ArrayList<Integer>();
		int tmp = 0;
		int p1;
		int p2;
		double score = 0.0d;
		double surfaceArea = 0.0d;
		
		voro.decomposite(data);
		faces = data.getFaces();
		amino = data.getAminos();
		
		System.out.println("\tdecomp done ... "+faces.size()+" faces");
		
		for(int id1: pepIds){
			if(faces.get(id1) == null){
				continue;
			}
			neighbors = faces.get(id1);
			surfaceFlag = false;
			surfaceArea = 0.0d;
			for(int id2 : neighbors.keySet()){
				if(gridIds.contains(id2) && neighbors.get(id2) > minContact){
					surfaceArea += neighbors.get(id2);
					surfaceFlag = true;
				}
			}
			if(surfaceFlag){
				//System.err.println(surfaceArea);
				surfaceIds.add(id1);
				if(surfaceArea > minContact){
					for(int id2 : neighbors.keySet()){
						if(!surfaceIds.contains(id2)){
							continue;
						}
						tmp = 0;
						for(int i = 25; i <= 150; i += 25){
							if(surfaceArea <= i*1.0d){
								break;
							}
							tmp++;
						}
						p1 = amino.get(id1).getOneLetterCode().charAt(0)-65;
						p2 = amino.get(id2).getOneLetterCode().charAt(0)-65;
						score += potential[p1][p2][tmp];
						
					}
				}
			}
		}
		
		return score;
	}
	
	/**
	 * TEST main method
	 */
	public static void main(String[] args){
		
		List<String> pdbIds = new ArrayList<String>();
		String pdbLoc = args[0];
		
		File file = new File(pdbLoc);
		for(File f : file.listFiles()){
			pdbIds.add(f.getName().substring(0, 7));
		}
		
		GridSolventPotential pot = new GridSolventPotential(args[1],pdbLoc, pdbIds, VoroPrepType.CA, Double.parseDouble(args[3]), Double.parseDouble(args[4]), Double.parseDouble(args[5]), Double.parseDouble(args[6]));
		pot.writeToFile(args[2]);
		System.out.println("done");
		
//		GridSolvensPotential pot = new GridSolvensPotential(args[2],args[1],VoroPrepType.CC);
//		PDBFileReader reader = new PDBFileReader(pdbLoc);
//		PDBEntry pdb = reader.readFromFolderById("1j2xA00");
//		System.out.println(pot.scoreModel(pdb));
		
	}

}
