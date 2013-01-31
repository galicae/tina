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

import bioinfo.energy.potential.preparation.voronoi.VoroPPWrap;
import bioinfo.energy.potential.preparation.voronoi.VoroPrepType;
import bioinfo.energy.potential.preparation.voronoi.VoroPrepare;
import bioinfo.energy.potential.preparation.voronoi.VoronoiData;
import bioinfo.proteins.AminoAcidName;
import bioinfo.proteins.DSSPEntry;
import bioinfo.proteins.DSSPFileReader;
import bioinfo.proteins.PDBEntry;
import bioinfo.proteins.PDBFileReader;

/**
 * just a dummy potential to evaluate the self-written solvensPotential
 * therefore scoreModel returns always 0.0d;
 * @author andreseitz
 *
 */
public class DSSPSolvensPotential implements IEnergy{

	/**
	 * potential contains the actual mean force potential
	 * [a][b][c]
	 * a contains aminoacid information (one letter code ascii - 65) of partner 1
	 * b contains aminoacid information (one letter code ascii - 65) of partner 2
	 * c contains area of face between the two partners with the following classes
	 * smaller then 25,50,75,100,125,150,bigger then 150, where all values smaller then 1 have to be ignored
	 */
	private double[][][] potential = new double[26][26][7];
	//private int[] aminoCount = new int[26];
	private final String pdbFolder;
	private final double MINCONTACT = 1.0d;
	private final int MINSURFACE = 10;
	private final double mkT = -0.582d;
	private final double gridHullExtend = 2.0d;
	private final double gridDensity = 1.0d;
	private final double gridClash = 4.0d;
	
	private final String[] mappingKeys = {"aminoacid1","aminoacid2","faceArea"};
	private String vorobin;
	private String tmpdir;
	
	public DSSPSolvensPotential(String vorobin, String dsspFolder, List<String> dsspIds){
		this.pdbFolder = dsspFolder;
		this.vorobin = vorobin;
		calculateFromDATA(dsspIds);
	}
	
	public DSSPSolvensPotential(String vorobin, String tmpdir, String dsspFolder, List<String> dsspIds){
		this.pdbFolder = dsspFolder;
		this.vorobin = vorobin;
		this.tmpdir = tmpdir;
		calculateFromDATA(dsspIds);

	}
	
	public DSSPSolvensPotential(String filename,String vorobin, String tmpdir){
		this.pdbFolder=null;
		this.vorobin = vorobin;
		this.tmpdir = tmpdir;
		this.readFromFile(filename);
	}
	
	public DSSPSolvensPotential(String filename,String vorobin){
		this.pdbFolder=null;
		this.vorobin = vorobin;
		this.readFromFile(filename);
	}
	
	@Override
	public void writeToFile(String filename) {
		BufferedWriter bw;
		try {
			bw = new BufferedWriter(new OutputStreamWriter(new FileOutputStream(filename)));
			
			for(int k = 0; k != potential[0][0].length; k++){
				bw.append("=="+k+"==\n");
				for(int i = 0; i != potential.length; i++){
					for(int j = 0; j != potential[0].length; j++){
						bw.append(potential[i][j][k]+"\t");
					}
					bw.append("\n");
				}
				bw.append("\n");
			}
			
			bw.flush();
			bw.close();
		} catch(Exception e){
			e.printStackTrace();
		}
	}

	@Override
	public void readFromFile(String filename) {
		try{
			BufferedReader br = new BufferedReader(new InputStreamReader(new FileInputStream(filename)));
			String line = null;
			Pattern classPattern = Pattern.compile("==(\\d)==");
			Matcher classMatcher = null;
			int k = 0;
			String[] fields = null;
			while((line = br.readLine())!= null){
				classMatcher = classPattern.matcher(line);
				if(classMatcher.find()){
					k = Integer.parseInt(classMatcher.group(1));
					for(int i = 0; i != 26; i++){
						line = br.readLine();
						fields = line.trim().split("\t");
						for(int j = 0; j != 26; j++){
							this.potential[i][j][k] = Double.parseDouble(fields[j]);
						}
					}
				}
			}
			
			br.close();
			
			
		}catch(Exception e){
			e.printStackTrace();
		}
	}

	@Override
	public void calculateFromDATA(List<String> dsspIds) {
		DSSPEntry dssp = null;
		DSSPFileReader reader = new DSSPFileReader(pdbFolder);
		VoronoiData data = null;
		VoroPrepare vprep = new VoroPrepare();
		Set<Integer> gridIds = null;
		Set<Integer> pepIds = null;
		HashMap<Integer,Integer> acc;
		HashMap<Integer, AminoAcidName> amino;
		HashMap<Integer, HashMap<Integer,Double>> faces;
		HashMap<Integer,Double> neighbors;
		HashMap<Integer, Double> neighborsNeighbors;
		boolean solvensNeighbor;
		boolean surfaceFlag = false;
		List<Integer> surfaceIds = new ArrayList<Integer>();
		int tmp = 0;
		int p1;
		int p2;
		int[] aminoCount = new int[26];
		int accessability = 0;
		
		for(String dsspId: dsspIds){
			
			System.out.println(">"+dsspId);
			
			//read from files, can be changed to db later!
			dssp = reader.readFromFolderById(dsspId);
			data = new VoronoiData(dsspId,null);
			data = vprep.reduceDSSP(dssp);
			acc = new HashMap<Integer,Integer>();
			for(int i = 0; i != dssp.getResIndex().length; i++){
				acc.put(dssp.getResIndex()[i], dssp.getAccesability()[i]);
			}
			pepIds = data.getAllInsertedIds();
			gridIds = data.fillGridWithoutClashes(data.extendHull(gridHullExtend), gridDensity, gridClash);
			VoroPPWrap voro;
			if(tmpdir == null){
				voro = new VoroPPWrap(vorobin);
			}else{
				voro = new VoroPPWrap(tmpdir,vorobin);
			}
			data = voro.decomposite(data);
			faces = data.getFaces();
			amino = data.getAmino();
			
			System.out.println("\tdecomp done ... "+faces.size()+" faces");
			
			for(int id1: pepIds){
				if(faces.get(id1) == null){
					continue;
				}
				neighbors = faces.get(id1);
				surfaceFlag = false;
				for(int id2 : neighbors.keySet()){
					if(gridIds.contains(id2) && neighbors.get(id2) > MINCONTACT){
						neighborsNeighbors = faces.get(id2);
						solvensNeighbor = false;
						for(int id3 : neighborsNeighbors.keySet()){
							if(gridIds.contains(id3) && neighborsNeighbors.get(id3) > MINCONTACT){
								solvensNeighbor = true;
							}
						}
						if(solvensNeighbor){
							surfaceFlag = true;
						}
					}
				}
				if(surfaceFlag){
					accessability = acc.get(id1);
					System.err.println(accessability);
					surfaceIds.add(id1);
					if(accessability > MINSURFACE){
						for(int id2 : neighbors.keySet()){
							if(pepIds.contains(id2) && neighbors.get(id2) > MINCONTACT){
								tmp = 0;
								for(int i = 25; i <= 150; i += 25){
									if(accessability <= i*1.0d){
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
		return 0.0d;
	}
	
	/**
	 * TEST main method
	 */
	public static void main(String[] args){
		
		List<String> dsspIds = new ArrayList<String>();
		String dsspLoc = args[0];
		
		File file = new File(dsspLoc);
		System.out.println(file.listFiles().length);
		for(File f : file.listFiles()){
			dsspIds.add(f.getName().substring(0, 7));
		}
		
		DSSPSolvensPotential pot = new DSSPSolvensPotential(args[1],dsspLoc, dsspIds);
		pot.writeToFile(args[2]);
		System.out.println("done");
		
//		GridSolvensPotential pot = new GridSolvensPotential(args[2],args[1],VoroPrepType.CC);
//		PDBFileReader reader = new PDBFileReader(pdbLoc);
//		PDBEntry pdb = reader.readFromFolderById("1j2xA00");
//		System.out.println(pot.scoreModel(pdb));
		
	}

}
