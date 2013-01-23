package bioinfo.energy.potential.preparation.voronoi;

import java.util.AbstractSet;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import bioinfo.energy.potential.preparation.voronoi.VoroPrepType;
import bioinfo.proteins.AminoAcidName;

public class VoronoiData {
	
	String id;
	VoroPrepType type;
	private HashMap<Integer, double[]> pointPosition;
	private HashMap<Integer, AminoAcidName> amino;
	private HashMap<Integer, HashMap<Integer,Double>> faces;
	private VoroDataState state;
	
	public VoronoiData(String id, VoroPrepType type){
		this.type = type;
		this.id = id;
		this.pointPosition = new HashMap<Integer, double[]>();
		this.amino = new HashMap<Integer,AminoAcidName>();
		this.faces = new HashMap<Integer,HashMap<Integer,Double>>();
		this.state = VoroDataState.INIT;
	}
	
	/**
	 * adds raw data to the object before calculation is started
	 * @param pos number of the residue
	 * @param coord coordinates the residue is placed at
	 * @param aa aminoacid name of the residue
	 */
	public void addValues(int pos, double[] coord, AminoAcidName aa){
		if(this.state == VoroDataState.INIT){
		pointPosition.put(pos,coord);
		amino.put(pos,aa);
		}else{
			return;
		}
	}
	
	/**
	 * adds data of voronoi decomposition to the object
	 * @param pos number of the residue
	 * @param neighbors array of numbers of neighboring residues
	 * @param areas size of areas shared with the neighbors named in neighbors
	 */
	public void addValues(int pos, int[] neighbors, double[] areas){
		HashMap<Integer,Double> area = new HashMap<Integer,Double>();
		for(int i = 0; i != neighbors.length; i++){
			area.put(neighbors[i], areas[i]);
		}
		this.faces.put(pos, area);
	}
	
	
	/**
	 * sets state of data from init to close
	 *no more data can be added to the object
	 */
	public void endInit(){
		if(this.state == VoroDataState.INIT){
			this.state = VoroDataState.CLOSE;
		}else{
			return;
		}
	}
	
	/**
	 * set state of data from close to calc
	 * which means calculation is done, values 
	 * are written ond nomare changes can be made
	 */
	public void endCalc(){
		if(this.state == VoroDataState.CLOSE){
			this.state = VoroDataState.CALC;
		}
	}
	
	/**
	 * 
	 * @return unsorted array of all coordinates
	 */
	public double[][] getAllCoord(){
		List<double[]> coord = new ArrayList<double[]>();
		for(int i : pointPosition.keySet()){
			coord.add(pointPosition.get(i));
		}
		return coord.toArray(new double[coord.size()][3]);
	}
	
	/**
	 * 
	 * @return String containing all information an format to pass the resulting file to voro++
	 */
	public String toVoroPPString(){
		String out = "";
		double[] coord;
		for(int i : pointPosition.keySet()){
			coord = pointPosition.get(i);
			out += i+"\t"+coord[0]+"\t"+coord[1]+"\t"+coord[2]+"\n";
		}
		return out;
		
	}
	
	public String getID(){
		return this.id;
	}
	
	public HashMap<Integer, HashMap<Integer,Double>> getFaces(){
		return faces;
	}
	
	public HashMap<Integer, AminoAcidName> getAmino(){
		return amino;
	}
	
	/**
	 * calculates min hull for voro++ voronoi decomposition, rounded to the next full numbers
	 * @return double[][] containing values: {{xmax,xmin},{ymax,ymin},{zmax,zmin}}
	 */
	public double[][] calcHull(){
		double[][] coord = this.getAllCoord(); //double[n][3]
		double[][] hull = {{1.0*Integer.MAX_VALUE, 1.0*Integer.MIN_VALUE},{1.0*Integer.MAX_VALUE, 1.0*Integer.MIN_VALUE},{1.0*Integer.MAX_VALUE, 1.0*Integer.MIN_VALUE}}; //x,z,y and min,max
		
		for(int n = 0; n != coord.length; n++){
			if(coord[n][0] < hull[0][0]){
				hull[0][0] = coord[n][0];
			}else if(coord[n][0] > hull[0][1]){
				hull[0][1] = coord[n][0];
			}
			
			if(coord[n][1] < hull[1][0]){
				hull[1][0] = coord[n][1];
			}else if(coord[n][1] > hull[1][1]){
				hull[1][1] = coord[n][1];
			}
			
			if(coord[n][2] < hull[2][0]){
				hull[2][0] = coord[n][2];
			}else if(coord[n][2] > hull[2][1]){
				hull[2][1] = coord[n][2];
			}
		}
		

		for(int i = 0; i != 3; i++){
			hull[i][0] = Math.floor(hull[i][0]);
			hull[i][1] = Math.ceil(hull[i][1]);
		}
		return hull;
	}
	
	/**
	 * 
	 * @param hull initial hull, calculated on points of peptide
	 * @param extension of hull eg 1.0 or 2.0 etc
	 * @return new hull extended by extension in each direction
	 */
	public double[][] extendHull (double extension){
		double[][] hull = this.calcHull();
		double[][] rhull = new double[3][2];
		for(int i = 0; i != 3; i++){
			rhull[i][0] = hull[i][0]-extension;
			rhull[i][1] = hull[i][1]+extension;
		}
		return rhull;
	}
	
	/**
	 * 
	 * @return Set of all existing identifiers in the dataset
	 */
	public Set<Integer> getAllInsertedIds(){
		return amino.keySet();
	}
	
	/**
	 * 
	 * @return next possible identifier, that makes sense given all existing identifiers
	 * it will be the max(given_identifiers)+1 value
	 */
	public int getNextIdentifier(){
		int max = 0;
		for(int i : amino.keySet()){
			if(i > max){
				max = i;
			}
		}
		return (max+1);
	}
	
	/**
	 * inserts lot of points in a grid like manner into the existing model to simulate water.
	 * clashing grid points will be removed
	 * @param hull which should be filled with water-like points
	 * @param gridDensity density of water-like points (size of water 2-2.5A, take 1.0 eg!)
	 * @param clashDistToNonGrid a water-like grid point will clash if euclidian distance is lower than given value
	 * normal distance between two points in non grid is about 3.8, so 2.0 should be a nice value to choose =)
	 * @return List of all grid ids 
	 */
	public Set<Integer> fillGridWithoutClashes(double[][] hull, double gridDensity, double clashDistToNonGrid){
		Set<Integer> gridIds = new HashSet<Integer>();
		Set<Integer> nonGridIds = pointPosition.keySet();
		int id = getNextIdentifier();
		double[] point;
		double[] tmp = null;
		boolean clash = false;
		for(double x = Math.floor(hull[0][0]); x <= Math.ceil(hull[0][1]); x = x+gridDensity){
			for(double y = Math.floor(hull[1][0]); y <= Math.ceil(hull[1][1]); y = y+gridDensity){
				for(double z = Math.floor(hull[2][0]); z <= Math.ceil(hull[2][1]); z = z+gridDensity){
					clash = false;
					for(int i : nonGridIds){
						point = pointPosition.get(i);
						if(Math.sqrt(((point[0]-x)*(point[0]-x))+((point[1]-y)*(point[1]-y))+((point[2]-z)*(point[2]-z))) < clashDistToNonGrid){
							clash = true;
							break;
						}
					}
					if(!clash){
						tmp = new double[3];
						tmp[0] = x;
						tmp[1] = y;
						tmp[2] = z;
						gridIds.add(id);
						addValues(id,tmp,AminoAcidName.U);
						id++;
					}
				}
			}
		}
		
		return gridIds;
	}
	
}
