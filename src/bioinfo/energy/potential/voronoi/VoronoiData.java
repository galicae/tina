package bioinfo.energy.potential.voronoi;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Set;

import bioinfo.energy.potential.voronoi.VoroPrepType;
import bioinfo.proteins.AminoAcid;
import bioinfo.proteins.AminoAcidName;
import bioinfo.proteins.Atom;
import bioinfo.proteins.AtomType;
import bioinfo.proteins.DSSPEntry;
import bioinfo.proteins.PDBEntry;

public class VoronoiData {

	private final String id;
	
	private HashMap<Integer, double[]> pointPositions;
	private HashMap<Integer, AminoAcidName> aminos;
	private HashMap<Integer, HashMap<Integer, Double>> faces;

	/**
	 * all peptide ids which were inserted
	 * subset of points keyset
	 */
	private Set<Integer> peptideIds;
	
	/**
	 * all grid ids which were inserted
	 * subset of points keyset
	 */
	private Set<Integer> gridIds;
	
	/**
	 * all grid ids which are directly accessable via other grid ids from the outside
	 * list is generated by method detectOuterGrid
	 */
	private Set<Integer> outerGridIds;

	
	
	
	/**
	 * simple constructor initialising all data
	 * @param id will be set as the data object id
	 */
	public VoronoiData(String id) {
		this.id = id;
		this.pointPositions = new HashMap<Integer, double[]>();
		this.aminos = new HashMap<Integer, AminoAcidName>();
		this.faces = new HashMap<Integer, HashMap<Integer, Double>>();

		this.peptideIds = new HashSet<Integer>();
		this.gridIds = new HashSet<Integer>();
		this.outerGridIds = new HashSet<Integer>();

	}

	/**
	 * adds raw data to the object before calculation is started
	 * 
	 * @param pos
	 *            number of the residue
	 * @param coord
	 *            coordinates the residue is placed at
	 * @param aa
	 *            aminoacid name of the residue
	 */
	public void addAminoValues(int pos, double[] coord, AminoAcidName aa) {
		pointPositions.put(pos, coord);
		aminos.put(pos, aa);
		peptideIds.add(pos);
	}

	/**
	 * adds raw data to the object before calculation is started
	 * 
	 * @param pos
	 *            number of the residue
	 * @param coord
	 *            coordinates the residue is placed at
	 * @param aa
	 *            aminoacid name of the residue
	 */
	public void addGridValues(int pos, double[] coord, AminoAcidName aa) {
		pointPositions.put(pos, coord);
		aminos.put(pos, aa);
		gridIds.add(pos);
	}

	/**
	 * adds data of voronoi decomposition to the object
	 * 
	 * @param pos
	 *            number of the residue
	 * @param neighbors
	 *            array of numbers of neighboring residues
	 * @param areas
	 *            size of areas shared with the neighbors named in neighbors
	 */
	public void addValues(int pos, int[] neighbors, double[] areas) {
		HashMap<Integer, Double> area = new HashMap<Integer, Double>();
		for (int i = 0; i != neighbors.length; i++) {
			area.put(neighbors[i], areas[i]);
		}
		this.faces.put(pos, area);
	}

	/**
	 * 
	 * @return unsorted array of all coordinates
	 */
	public double[][] getAllCoord() {
		List<double[]> coord = new ArrayList<double[]>();
		for (int i : pointPositions.keySet()) {
			coord.add(pointPositions.get(i));
		}
		return coord.toArray(new double[coord.size()][3]);
	}

	/**
	 * 
	 * @return String containing all information an format to pass the resulting
	 *         file to voro++
	 */
	public String toVoroPPString() {
		String out = "";
		double[] coord;
		for (int i : pointPositions.keySet()) {
			coord = pointPositions.get(i);
			out += i + "\t" + coord[0] + "\t" + coord[1] + "\t" + coord[2] + "\n";
		}
		return out;

	}

	/**
	 * 
	 * @return the ID of the voronoi data set or the dssp/pdb model 
	 */
	public String getID() {
		return this.id;
	}

	/**
	 * 
	 * @return the face-areas calculated by Voronoi Decomposition
	 */
	public HashMap<Integer, HashMap<Integer, Double>> getFaces() {
		return faces;
	}

	/**
	 * 
	 * @return all inserted points
	 */
	public HashMap<Integer, AminoAcidName> getAminos() {
		return aminos;
	}

	/**
	 * 
	 * @return all peptide Ids inserted into the data object
	 */
	public Set<Integer> getPepIds() {
		return this.peptideIds;
	}

	/**
	 * 
	 * @return all grid ids inserted into the data object
	 */
	public Set<Integer> getGridIds() {
		return this.gridIds;
	}
	
	/**
	 * @return all outer grid ids which means grid ids which are accessable via all other outer grid ids
	 */
	public Set<Integer> getOuterGridIds(){
		return outerGridIds;
	}
	
	/**
	 * calculates min hull for voro++ voronoi decomposition, rounded to the next
	 * full numbers
	 * 
	 * @return double[][] containing values:
	 *         {{xmax,xmin},{ymax,ymin},{zmax,zmin}}
	 */
	public double[][] calcHull() {
		double[][] coord = this.getAllCoord(); // double[n][3]
		//X,Y,Z and MIN,MAX
		double[][] hull = { { 1.0 * Integer.MAX_VALUE, 1.0 * Integer.MIN_VALUE }, { 1.0 * Integer.MAX_VALUE, 1.0 * Integer.MIN_VALUE }, { 1.0 * Integer.MAX_VALUE, 1.0 * Integer.MIN_VALUE } };
		
		for (int n = 0; n != coord.length; n++) {
			if (coord[n][0] < hull[0][0]) {
				hull[0][0] = coord[n][0];
			} else if (coord[n][0] > hull[0][1]) {
				hull[0][1] = coord[n][0];
			}

			if (coord[n][1] < hull[1][0]) {
				hull[1][0] = coord[n][1];
			} else if (coord[n][1] > hull[1][1]) {
				hull[1][1] = coord[n][1];
			}

			if (coord[n][2] < hull[2][0]) {
				hull[2][0] = coord[n][2];
			} else if (coord[n][2] > hull[2][1]) {
				hull[2][1] = coord[n][2];
			}
		}

		for (int i = 0; i != 3; i++) {
			hull[i][0] = Math.floor(hull[i][0]);
			hull[i][1] = Math.ceil(hull[i][1]);
		}
		return hull;
	}

	/**
	 * 
	 * @param hull
	 *            initial hull, calculated on points of peptide
	 * @param extension
	 *            of hull eg 1.0 or 2.0 etc
	 * @return new hull extended by extension in each direction
	 */
	private double[][] extendHull(double extension) {
		double[][] hull = this.calcHull();
		double[][] rhull = new double[3][2];
		for (int i = 0; i != 3; i++) {
			rhull[i][0] = hull[i][0] - extension;
			rhull[i][1] = hull[i][1] + extension;
		}
		return rhull;
	}

	/**
	 * 
	 * @return next possible identifier, that makes sense given all existing
	 *         identifiers it will be the max(given_identifiers)+1 value
	 */
	public int getNextIdentifier() {
		int max = 0;
		for (int i : aminos.keySet()) {
			if (i > max) {
				max = i;
			}
		}
		return (max + 1);
	}

	/**
	 * inserts lot of points in a grid like manner into the existing model to
	 * simulate water. clashing grid points will be removed
	 * 
	 * @param hull
	 *            which should be filled with water-like points
	 * @param gridDensity
	 *            density of water-like points (size of water 2-2.5A, take 1.0
	 *            eg!)
	 * @param clashDistToNonGrid
	 *            a water-like grid point will clash if euclidian distance is
	 *            lower than given value normal distance between two points in
	 *            non grid is about 3.8, so 2.0 should be a nice value to choose
	 *            =)
	 */
	public void fillGridWithoutClashes(double gridExtend, double gridDensity, double clashDistToNonGrid) {
		int id = getNextIdentifier();
		double[] point;
		double[] tmp = null;
		boolean clash = false;
		double[][] hull = this.extendHull(gridExtend);
		for (double x = Math.floor(hull[0][0]); x <= Math.ceil(hull[0][1]); x = x + gridDensity) {
			for (double y = Math.floor(hull[1][0]); y <= Math.ceil(hull[1][1]); y = y + gridDensity) {
				for (double z = Math.floor(hull[2][0]); z <= Math.ceil(hull[2][1]); z = z + gridDensity) {
					clash = false;
					for (int i : pointPositions.keySet()) {
						point = pointPositions.get(i);
						if (Math.sqrt(((point[0] - x) * (point[0] - x)) + ((point[1] - y) * (point[1] - y)) + ((point[2] - z) * (point[2] - z))) < clashDistToNonGrid) {
							clash = true;
							break;
						}
					}
					if (!clash) {
						tmp = new double[3];
						tmp[0] = x;
						tmp[1] = y;
						tmp[2] = z;
						addGridValues(id, tmp, AminoAcidName.U);
						id++;
					}
				}
			}
		}
	}

	/**
	 * reduces pdb entry to voronoi data object
	 * @param type VoroPrpType to be used
	 * @param pdb entry that will be represented by the data object
	 */
	public void reducePDB(VoroPrepType type, PDBEntry pdb) {
		switch (type) {
		case CA: {
			Atom tmp = null;
			for (int i = 0; i != pdb.length(); i++) {
				if ((tmp = pdb.getAminoAcid(i).getAtomByType(AtomType.CA)) != null) {
					addAminoValues(i, tmp.getPosition(), pdb.getAminoAcid(i).getName());
				}
			}
		}
			;
			break;
		case CC: {
			AminoAcid tmp;
			Atom temp;
			double[][] coord;
			double[] centroid;
			for (int i = 0; i != pdb.length(); i++) {
				tmp = pdb.getAminoAcid(i);
				coord = new double[tmp.getAtomNumber()][3];
				for (int j = 0; j != tmp.getAtomNumber(); j++) {
					temp = tmp.getAtom(j);
					coord[j] = temp.getPosition();
				}
				centroid = calculateCentroid(coord);
				addAminoValues(i, centroid, tmp.getName());
			}
		}
			;
			break;
		case SC: {
			AminoAcid aaTmp;
			Atom aTmp;
			List<double[]> coordTmp;
			double[] centroid;
			for (int i = 0; i != pdb.length(); i++) {
				aaTmp = pdb.getAminoAcid(i);
				coordTmp = new ArrayList<double[]>();
				for (int j = 0; j != aaTmp.getAtomNumber(); j++) {
					aTmp = aaTmp.getAtom(j);
					if (aTmp.getType() != AtomType.C && aTmp.getType() != AtomType.CA && aTmp.getType() != AtomType.O && aTmp.getType() != AtomType.N) {
						coordTmp.add(aTmp.getPosition());
					}
				}
				centroid = calculateCentroid(coordTmp.toArray(new double[coordTmp.size()][3]));
				addAminoValues(i, centroid, aaTmp.getName());
			}
		}
			break;
		}
	}

	/**
	 * reduces dssp entry to voronoi data object
	 * @param dssp dssp entry tht will be represented by the data object
	 */
	public void reduceDSSP(DSSPEntry dssp) {

		for (int i = 0; i != dssp.getLength(); i++) {
			addAminoValues(i, dssp.getCaTrace()[i], dssp.getNames()[i]);
		}

	}

	/**
	 * just calculate the centroid of some points
	 * @param coord[][] with nx3 containing all coordinates of points whose centroid should be calculated
	 * @return the centroid coordinates of all inserted points
	 */
	private double[] calculateCentroid(double[][] coord) {
		double cx = 0, cy = 0, cz = 0;
		for (int i = 0; i < coord.length; i++) {
			cx += coord[i][0];
			cy += coord[i][1];
			cz += coord[i][2];
		}
		cx = cx / (coord.length * 1.0);
		cy = cy / (coord.length * 1.0);
		cz = cz / (coord.length * 1.0);
		double[] result = { cx, cy, cz };
		return result;
	}

	/**
	 * detects outer grid points and writes them to outerGridIds
	 */
	public void detectOuterGrid(double minContact){
		List<Integer> notVisited = new LinkedList<Integer>();
		boolean isPartOfBox = false;
		for(int id : gridIds){
			isPartOfBox = false;
			if(faces.get(id) != null){
				for(int neighbor : faces.get(id).keySet()){
					if(neighbor < 0){
						isPartOfBox = true;
						break;
					}
				}
				if(isPartOfBox){
					notVisited.add(id);
					//outerGridIds.add(id);
				}
			}
		}
		
		int current;
		while(!notVisited.isEmpty()){
			current = notVisited.get(0);
			notVisited.remove(0);
			if(!outerGridIds.contains(current) && gridIds.contains(current)){
				outerGridIds.add(current);
				for(int neighbor : faces.get(current).keySet()){
					if(faces.get(current).get(neighbor) > minContact){
						notVisited.add(neighbor);
					}
				}
			}
		}
	}
	
}