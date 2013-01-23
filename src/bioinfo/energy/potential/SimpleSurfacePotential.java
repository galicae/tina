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
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import bioinfo.energy.potential.preparation.voronoi.VoroPPWrap;
import bioinfo.energy.potential.preparation.voronoi.VoroPrepType;
import bioinfo.energy.potential.preparation.voronoi.VoroPrepare;
import bioinfo.energy.potential.preparation.voronoi.VoronoiData;
import bioinfo.proteins.AminoAcidName;
import bioinfo.proteins.PDBEntry;
import bioinfo.proteins.PDBFileReader;

public class SimpleSurfacePotential implements IEnergy{
	
	/**
	 * potential contains the actual mean force potential
	 * [a][b][c]
	 * a contains aminoacid information (one letter code ascii - 65) of partner 1
	 * b contains aminoacid information (one letter code ascii - 65) of partner 2
	 * c contains area of face between the two partners with the following classes
	 * smaller then 4,8,16,32,64,128,256,512,bigger then 512, where all values smaller then 2 have to be ignored
	 */
	private double[][][] potential = new double[26][26][9];
	private final String pdbFolder;
	private final VoroPrepType type;
	private final double MINCONTACT = 2.0d;
	private final double mkT = -0.582d;
	private final String[] mappingKeys = {"aminoacid1","aminoacid2","faceArea"};
	private String vorobin = null;
	private String tmpdir = null;
	
	public SimpleSurfacePotential(String vorobin, String pdbFolder, List<String> pdbIds, VoroPrepType type){
		this.pdbFolder = pdbFolder;
		this.type = type;
		calculateFromPDBs(pdbIds);
		this.vorobin = vorobin;
	}
	
	public SimpleSurfacePotential(String vorobin, String tmpdir, String pdbFolder, List<String> pdbIds, VoroPrepType type){
		this.pdbFolder = pdbFolder;
		this.type = type;
		calculateFromPDBs(pdbIds);
		this.vorobin = vorobin;
		this.tmpdir = tmpdir;
	}
	
	public SimpleSurfacePotential(String filename){
		this.type=null;
		this.pdbFolder=null;
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
	public void calculateFromPDBs(List<String> pdbIds) {
		PDBEntry pdb = null;
		PDBFileReader reader = new PDBFileReader(pdbFolder);
		VoronoiData data = null;
		VoroPrepare vprep = new VoroPrepare();
		HashMap<Integer, AminoAcidName> amino;
		HashMap<Integer, HashMap<Integer,Double>> faces;
		HashMap<Integer,Double> neighbors;
		boolean surfaceFlag = false;
		List<Integer> surfaceIds = new ArrayList<Integer>();
		int tmp = 0;
		int p1;
		int p2;
		int count = 0;
		
		
		for(String pdbId: pdbIds){
			
			System.out.println(">"+pdbId);
			
			//read from files, can be changed to db later!
			pdb = reader.readFromFolderById(pdbId);
			data = new VoronoiData(pdbId, type);
			data = vprep.reducePDB(type, pdb);
			VoroPPWrap voro;
			if(this.tmpdir == null){
				voro = new VoroPPWrap(vorobin);
			}else{
				voro = new VoroPPWrap(vorobin,tmpdir);
			}
			data = voro.decomposite(data);
			faces = data.getFaces();
			amino = data.getAmino();
			
			System.out.println("\tdecomp done ... "+faces.size()+" faces");
			
			for(int id1: amino.keySet()){
				if(faces.get(id1) == null){
					continue;
				}
				neighbors = faces.get(id1);
				surfaceFlag = false;
				for(int id2 : neighbors.keySet()){
					if(id2 < 0 && neighbors.get(id2) > MINCONTACT){
						surfaceFlag = true;
						surfaceIds.add(id1);
						break;
					}
				}
				if(surfaceFlag){
					for(int id2 : neighbors.keySet()){
						if(id2 > 0 && neighbors.get(id2) > MINCONTACT){
							tmp = 0;
							for(int i = 4; i <= 512; i *= 2){
								if(neighbors.get(id2) <= i*1.0d){
									break;
								}
								tmp++;
							}
							p1 = amino.get(id1).getOneLetterCode().charAt(0)-65;
							p2 = amino.get(id2).getOneLetterCode().charAt(0)-65;
							potential[p1][p2][tmp]++;
							potential[p2][p1][tmp]++;
							count++;
						}
					}
				}
			}
			
		}
		
		for(int i = 0; i != 26; i++){
			for(int j = 0; j != 26; j++){
				for(int k = 0; k != 9; k++){
					potential[i][j][k] = mkT*Math.log(potential[i][j][k]+1/count+1);
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
		for(int i = 4; i <= 512; i *= 2){
			if(area <= i*1.0d){
				break;
			}
			k++;
		}
		
		int i = aa1-65;
		int j = aa2-65;
		
		return potential[i][j][k];
		
	}
	
	
	
	
	
	/**
	 * TEST main method
	 */
	public static void main(String[] args){
		
		SimpleSurfacePotential pot1 = new SimpleSurfacePotential("/Users/andreseitz/Desktop/surfacePotential_CC.pot");
		System.exit(1);
		List<String> pdbIds = new ArrayList<String>();
		String pdbLoc = "/Users/andreseitz/Documents/uni/CIP/gobi/STRUCTURES/";
		File file = new File(pdbLoc);
		for(File f : file.listFiles()){
			pdbIds.add(f.getName().substring(0, 7));
		}
		SimpleSurfacePotential pot = new SimpleSurfacePotential("./tools/voro++",pdbLoc, pdbIds, VoroPrepType.CC);
		pot.writeToFile("/Users/andreseitz/Desktop/surfacePotential_CC.pot");
		System.out.println("done");
		
	}

	@Override
	public String[] getMapKeys() {
		return this.mappingKeys;
	}

}
