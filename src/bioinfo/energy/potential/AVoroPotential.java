package bioinfo.energy.potential;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.InputStreamReader;
import java.io.OutputStreamWriter;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import bioinfo.energy.potential.preparation.voronoi.VoroPPWrap;
import bioinfo.energy.potential.preparation.voronoi.VoroPrepType;
import bioinfo.energy.potential.preparation.voronoi.VoronoiData;
import bioinfo.proteins.DSSPEntry;
import bioinfo.proteins.PDBEntry;
import bioinfo.proteins.PDBFileReader;

public abstract class AVoroPotential implements IEnergy{

	protected double[][][] potential;
	protected String[] mappingKeys;
	
	protected VoronoiData data;
	protected VoroPPWrap voro;
	
	public AVoroPotential(String vorobin, String tmpdir){
		voro = new VoroPPWrap(tmpdir,vorobin);
	}
	
	public AVoroPotential(String vorobin){
		voro = new VoroPPWrap("/tmp/", vorobin);
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
	
	public VoronoiData prepareSimple(PDBEntry pdb, VoroPrepType type){
		data = new VoronoiData(pdb.getID());
		data.reducePDB(type, pdb);
		voro.decomposite(data);
		return data;
	}
	
	public VoronoiData prepareSimple(DSSPEntry dssp){
		data = new VoronoiData(dssp.getID());
		data.reduceDSSP(dssp);
		voro.decomposite(data);
		return data;
	}
	
	public VoronoiData prepareWithGrid(PDBEntry pdb, VoroPrepType type, double gridExtend, double gridDensity, double gridClash, double minContact){
		data = new VoronoiData(pdb.getID());
		data.reducePDB(type, pdb);
		data.fillGridWithoutClashes(gridExtend, gridDensity, gridClash);
		voro.decomposite(data);
		data.detectOuterGrid(minContact);
		return data;
	}
	
	public VoronoiData prepareWithGrid(DSSPEntry dssp, double gridExtend, double gridDensity, double gridClash, double minContact){
		data = new VoronoiData(dssp.getID());
		data.reduceDSSP(dssp);
		data.fillGridWithoutClashes(gridExtend, gridDensity, gridClash);
		voro.decomposite(data);
		data.detectOuterGrid(minContact);
		return data;
	}	
	

}
