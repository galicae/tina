package bioinfo.energy.potential;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
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
import bioinfo.energy.potential.preparation.voronoi.VoronoiData;
import bioinfo.proteins.AminoAcidName;
import bioinfo.proteins.DSSPEntry;
import bioinfo.proteins.DSSPFileReader;
import bioinfo.proteins.PDBEntry;

/**
 * just a dummy potential to evaluate the self-written solvensPotential
 * therefore scoreModel returns always 0.0d;
 * @author andreseitz
 *
 */
public class DSSPSolventPotential extends AVoroPotential{

	/**
	 * potential contains the actual mean force potential
	 * [a][b][c]
	 * a contains aminoacid information (one letter code ascii - 65) of partner 1
	 * b contains aminoacid information (one letter code ascii - 65) of partner 2
	 * c contains area of face between the two partners with the following classes
	 * smaller then 25,50,75,100,125,150,bigger then 150, where all values smaller then 1 have to be ignored
	 */
	private final String dsspFolder;
	private final double minContact;
	private final double mkT = -0.582d;
	private final double gridHullExtend;
	private final double gridDensity;
	private final double gridClash;
	
	private final String[] mappingKeys = {"aminoacid1","aminoacid2","faceArea"};
	private String vorobin;
	private String tmpdir;
	
	public DSSPSolventPotential(String vorobin, String dsspFolder, List<String> dsspIds){
		super(vorobin);
		this.potential = new double[26][26][7];
		this.dsspFolder = dsspFolder;
		
		this.minContact = 1.0d;
		this.gridHullExtend = 2.0d;
		this.gridDensity = 3.0d;
		this.gridClash = 4.0d;
		
		calculateFromDATA(dsspIds);
	}
	
	public DSSPSolventPotential(String vorobin, String tmpdir, String dsspFolder, List<String> dsspIds){
		super(vorobin,tmpdir);
		this.potential = new double[26][26][7];
		this.dsspFolder = dsspFolder;
		
		this.minContact = 1.0d;
		this.gridHullExtend = 2.0d;
		this.gridDensity = 3.0d;
		this.gridClash = 4.0d;
		
		calculateFromDATA(dsspIds);

	}
	
	public DSSPSolventPotential(String vorobin, String dsspFolder, List<String> dsspIds, double minContact, double gridHullExtend, double gridDensity, double gridClash){
		super(vorobin);
		this.potential = new double[26][26][7];
		this.dsspFolder = dsspFolder;
		
		this.minContact = minContact;
		this.gridHullExtend = gridHullExtend;
		this.gridDensity = gridDensity;
		this.gridClash = gridClash;
		
		calculateFromDATA(dsspIds);
	}
	
	public DSSPSolventPotential(String vorobin, String tmpdir, String dsspFolder, List<String> dsspIds, double minContact, double gridHullExtend, double gridDensity, double gridClash){
		super(vorobin, tmpdir);
		this.potential = new double[26][26][7];
		this.dsspFolder = dsspFolder;
		
		this.minContact = minContact;
		this.gridHullExtend = gridHullExtend;
		this.gridDensity = gridDensity;
		this.gridClash = gridClash;
		
		calculateFromDATA(dsspIds);

	}
	
	public DSSPSolventPotential(String filename,String vorobin, String tmpdir){
		super(vorobin, tmpdir);
		this.potential = new double[26][26][7];
		this.dsspFolder = null;
		
		this.minContact = 1.0d;
		this.gridHullExtend = 2.0d;
		this.gridDensity = 3.0d;
		this.gridClash = 4.0d;
		
		this.readFromFile(filename);
	}
	
	public DSSPSolventPotential(String filename,String vorobin){
		super(vorobin);
		this.potential = new double[26][26][7];
		this.dsspFolder = null;
		
		this.minContact = 1.0d;
		this.gridHullExtend = 2.0d;
		this.gridDensity = 3.0d;
		this.gridClash = 4.0d;
		
		this.readFromFile(filename);
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
		HashMap<Integer, HashMap<Integer,Double>> faces;
		HashMap<Integer,Double> neighbors;
		boolean surfaceFlag = false;
		List<Integer> surfaceIds = new ArrayList<Integer>();
		int tmp = 0;
		int p1;
		int p2;
		int[] aminoCount = new int[26];
		int accessability = 0;
		double surfaceArea = 0.0d;
		
		for(String dsspId: dsspIds){
			//test->
			BufferedWriter bw = null;
			try {
				bw = new BufferedWriter(new OutputStreamWriter(new FileOutputStream("/Users/andreseitz/Desktop/"+dsspId+".res")));
			} catch (FileNotFoundException e) {
				e.printStackTrace();
			}
			//->test
			
			
			System.out.println(">"+dsspId);
			
			//read from files, can be changed to db later!
			dssp = reader.readFromFolderById(dsspId);
			
			//test->
			AminoAcidName[] names = dssp.getNames();
			int[] index = dssp.getResIndex();
			HashMap<Integer, Double> surfaces = new HashMap<Integer,Double>();
			//->test
			
		
			data = this.prepareWithGrid(dssp, gridHullExtend, gridDensity, gridClash, minContact);
			acc = dssp.getAccesability();
			pepIds = data.getPepIds();
			solventIds = data.getOuterGridIds();
			faces = data.getFaces();
			point = data.getPoints();
			
			
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
					accessability = acc[id1];
					//test->
					surfaces.put(id1,surfaceArea);
					//->test
					System.err.println(surfaceArea+"\t"+accessability);
					surfaceIds.add(id1);
					if(accessability > minContact){
						for(int id2 : neighbors.keySet()){
							if(!surfaceIds.contains(id2)){
								continue;
							}
							tmp = 0;
							for(int i = 25; i <= 150; i += 25){
								if(accessability <= i*1.0d){
									break;
								}
								tmp++;
							}
			
							p1 = point.get(id1).getOneLetterCode().charAt(0)-65;
							p2 = point.get(id2).getOneLetterCode().charAt(0)-65;
							potential[p1][p2][tmp]++;
							potential[p2][p1][tmp]++;
							aminoCount[p1]++;
						}
					}
				}
			}
			
			//test->
			for(int i = 0; i != pepIds.size(); i++){
				try{
					if(surfaces.containsKey(i)){
						bw.append(names[i]+"\t"+index[i]+"\t"+acc[i]+"\t"+surfaces.get(i)+"\t"+(surfaces.get(i)-acc[i])+"\n");
					}else{
						bw.append(names[i]+"\t"+index[i]+"\t"+acc[i]+"\t"+0+"\t"+((-1)*acc[i])+"\n");	
					}
				}catch(Exception e){
					e.printStackTrace();
				}
			}
			try{
				bw.flush();
				bw.close();
			}catch(Exception e){
				e.printStackTrace();
			}
			//->test
			
		}
		
		for(int i = 0; i != 26; i++){
			for(int j = 0; j != 26; j++){
				for(int k = 0; k != 7; k++){
					if(i == j){
						potential[i][j][k] = mkT*Math.log(((potential[i][j][k]+1))/(((aminoCount[i]+aminoCount[j])/2)+1));
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
		String vorobin = args[1];
		String outfile = args[2];
		String minContact = args[3];
		String gridExtend = args[4];
		String gridDensity = args[5];
		String gridClash = args[6];
		
		File file = new File(dsspLoc);
		for(File f : file.listFiles()){
			dsspIds.add(f.getName().substring(0, 7));
		}
		
		DSSPSolventPotential pot = new DSSPSolventPotential(vorobin, dsspLoc, dsspIds ,Double.parseDouble(minContact), Double.parseDouble(gridExtend), Double.parseDouble(gridDensity), Double.parseDouble(gridClash));
		pot.writeToFile(outfile);
		System.out.println("done");
		
//		GridSolvensPotential pot = new GridSolvensPotential(args[2],args[1],VoroPrepType.CC);
//		PDBFileReader reader = new PDBFileReader(pdbLoc);
//		PDBEntry pdb = reader.readFromFolderById("1j2xA00");
//		System.out.println(pot.scoreModel(pdb));
		
	}

}
