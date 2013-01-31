package db.mysql;

import java.sql.PreparedStatement;
import java.sql.ResultSet;
import java.sql.SQLException;
import java.sql.Statement;
import java.util.ArrayList;
import java.util.List;

import bioinfo.proteins.AminoAcid;
import bioinfo.proteins.AminoAcidName;
import bioinfo.proteins.Atom;
import bioinfo.proteins.DSSPEntry;
import bioinfo.proteins.PDBEntry;
import bioinfo.proteins.SecStructEight;

public class DBConnector extends MysqlWrapper{

	private static final String tablename = "pdb";
	
	private static final String[] pdbfields = {"id","pdb_id","chainID", "chainIDNum","length"};
	private static final String[] dsspfields = {"id","ss","sa", "phi","psi","x_ca","y_ca","z_ca","aminoacid_id","pdb_id"};
	private static final String[] aafields = {"id","name","res_index","numberofAtom","pdb_id"};
	private static final String[] atomfields = {"id","type","x","y","z","aminoacid_id"};
	
	private static final String setPDBEntry = "insert into "+tablename+" ("+pdbfields[1]+","+pdbfields[2]+","+pdbfields[3]+","+pdbfields[4]+") values (?,?,?,?)";	
	private static final String getPDBById = "select * from pdb join aminoacid on pdb.id = aminoacid.pdb_id join atom on aminoacid.id = atom.aminoacid_id" +
			" where pdb.pdb_id = ? and pdb.chainID = ? and pdb.chainIDNum = ?";
	
	private static final String setDSSPEntry = "insert into dssp ("
			+ dsspfields[1] + "," + dsspfields[2] + "," + dsspfields[3] + ","
			+ dsspfields[4] + "," + dsspfields[5] + "," + dsspfields[6] + ","
			+ dsspfields[7] + "," + dsspfields[8] + "," + dsspfields[9] + ") values (?,?,?,?,?,?,?,?,?)";	
	private static final String getDSSPById = "select dssp.*,aminoacid.name,aminoacid.res_index from dssp join pdb on pdb.id = dssp.pdb_id join aminoacid on aminoacid.id = dssp.aminoacid_id" +
			" where pdb.pdb_id = ? and pdb.chainID = ? and pdb.chainIDNum = ?";
	
	public DBConnector(MysqlDBConnection connection) {
		super(connection);
	}
	
	@Override
	String getTablename() {
		return tablename;
	}

	@Override
	String[] getFields() {
		return pdbfields;
	}
	
	public PDBEntry getPDB(String id){
		PreparedStatement stmt = connection.createStatement(getPDBById);
		List<Atom> atoms = new ArrayList<Atom>();
		List<AminoAcid> aminos = new ArrayList<AminoAcid>();
		
		String chainID = id.substring(4,5);
		int chainIDNum=Integer.parseInt(id.substring(5));
		String pdbid = id.substring(0, 4);
		
		try{
			stmt.setString(1, pdbid);
			stmt.setString(2, chainID);
			stmt.setInt(3, chainIDNum);
			ResultSet res = stmt.executeQuery();
			
			int atomquantity;
			double[] pos_temp;
			
			if(pdbExist(id) != -1){
				while(res.next()){
					//read out aminos and corresponding atoms
					atomquantity = res.getInt(aafields[3]);
					pos_temp = new double[3];
					pos_temp[0] = res.getDouble(atomfields[2]);
					pos_temp[1] = res.getDouble(atomfields[3]);
					pos_temp[2] = res.getDouble(atomfields[4]);
					atoms.add(new Atom(res.getString(atomfields[1]),pos_temp));
					
					for (int j = 0; j < atomquantity-1; j++) {
						res.next();
						pos_temp = new double[3];
						pos_temp[0] = res.getDouble(atomfields[2]);
						pos_temp[1] = res.getDouble(atomfields[3]);
						pos_temp[2] = res.getDouble(atomfields[4]);
						atoms.add(new Atom(res.getString(atomfields[1]),pos_temp));
					}
					aminos.add(new AminoAcid(AminoAcidName.getAAFromOLC(res.getString(aafields[1])),res.getInt(aafields[2]),atoms.toArray(new Atom[atomquantity])));
					atoms.clear();
				}
				return new PDBEntry(id,aminos.toArray(new AminoAcid[aminos.size()]));
			} else {
				return null;
			}
		}catch(SQLException e){
			e.printStackTrace();
			return null;
		}
	}
	
	public DSSPEntry getDSSP(String id){
		PreparedStatement stmt = connection.createStatement(getDSSPById);
		List<AminoAcidName> aminos = new ArrayList<AminoAcidName>();
		List<Integer> resIndex = new ArrayList<Integer>();
		List<SecStructEight> secStruct = new ArrayList<SecStructEight>();
		List<double[]> caTrace = new ArrayList<double[]>();
		List<Integer> sa = new ArrayList<Integer>();
		List<Double> phi = new ArrayList<Double>();
		List<Double> psi = new ArrayList<Double>();
		
		String pdbid = id.substring(0, 4);
		String chainID = id.substring(4,5);
		int chainIDNum=Integer.parseInt(id.substring(5));
		
		try{
			stmt.setString(1, pdbid);
			stmt.setString(2, chainID);
			stmt.setInt(3, chainIDNum);
			ResultSet res = stmt.executeQuery();

			if(pdbExist(id) != -1){		
				while(res.next()){
					//read out dssp
					aminos.add(AminoAcidName.getAAFromOLC(res.getString("aminoacid."+aafields[1])));
					resIndex.add(res.getInt("aminoacid."+aafields[2]));
					secStruct.add(SecStructEight.getSSFromChar(res.getString(dsspfields[1]).charAt(0)));
					caTrace.add(new double[] {res.getDouble(dsspfields[5]),res.getDouble(dsspfields[6]),res.getDouble(dsspfields[7])});
					sa.add(res.getInt(dsspfields[2]));
					phi.add(res.getDouble(dsspfields[3]));
					psi.add(res.getDouble(dsspfields[4]));
					
				}
				return new DSSPEntry(id,aminos,resIndex,secStruct,sa,phi,psi,caTrace);
			} else {
				return null;
			}
		}catch(SQLException e){
			e.printStackTrace();
			return null;
		}
	}
	
	public int pdbExist(String id){
		String query = "select id from "+tablename+" where "+pdbfields[1]+" = ? and "+pdbfields[2]+" = ? and "+pdbfields[3]+" = ?";
		PreparedStatement stmt = connection.createStatement(query);
		
		String pdbid = id.substring(0,4);
		String chain = id.substring(4,5);
		int chainIdNum = Integer.parseInt(id.substring(5));
		try{
			stmt.setString(1, pdbid);
			stmt.setString(2, chain);
			stmt.setInt(3, chainIdNum);
			
			ResultSet res = stmt.executeQuery();
			
			if(res.first()){
				return res.getInt(pdbfields[0]);
			}
			return -1;
		}catch(SQLException e){
			e.printStackTrace();
			return -1;
		}
	}
	
	public int dsspExist(String id){
		String query = "select dssp.id from dssp join pdb on pdb.id = dssp.pdb_id where pdb."+pdbfields[1]+" = ? and "+pdbfields[2]+" = ? and "+pdbfields[3]+" = ?";
		PreparedStatement stmt = connection.createStatement(query);
		
		String pdbid = id.substring(0,4);
		String chain = id.substring(4,5);
		int chainIdNum = Integer.parseInt(id.substring(5));
		try{
			stmt.setString(1, pdbid);
			stmt.setString(2, chain);
			stmt.setInt(3, chainIdNum);
			
			ResultSet res = stmt.executeQuery();
			
			if(res.first()){
				return res.getInt(pdbfields[0]);
			}
			return -1;
		}catch(SQLException e){
			e.printStackTrace();
			return -1;
		}
	}
	
	public String[] getIDs(){
		Statement stmt = connection.createStatement();
		try {
			ResultSet res = stmt.executeQuery("select id from "+getTablename());
			List<String> ids = new ArrayList<String>();
			while(res.next()){
				ids.add(res.getString(pdbfields[0]));
			}
			return ids.toArray(new String[ids.size()]);
		} catch (SQLException e) {
			e.printStackTrace();
			return null;
		}
	}
		
	private int getLastId(){
		int lastid;
		Statement stmt = connection.createStatement();
		ResultSet res;
		try {
			res = stmt.executeQuery("Select LAST_INSERT_ID() from "+tablename);
			res.first();
			lastid = res.getInt(1);
		} catch (SQLException e) {
			e.printStackTrace();
			return 0;
		}	
		return lastid;
	}
	
	public boolean addPDBEntry(PDBEntry entry){		
		
		if(pdbExist(entry.getID()+entry.getChainID()+entry.getChainIDNum()) == -1){
			try{
				AAConnector aaconnector = new AAConnector(this.connection);
				PreparedStatement stmt = connection.createStatement(setPDBEntry);
				AminoAcid amino;
	
				//insert pdbentry
				stmt.setString(1,entry.getID());
				stmt.setString(2,String.valueOf(entry.getChainID()));
				stmt.setInt(3,entry.getChainIDNum());
				stmt.setInt(4, entry.length());
				stmt.execute();
				int pdbid = getLastId();
				
				for (int i = 0; i < entry.length(); i++) {
					//insert amino
					amino = entry.getAminoAcid(i);
					aaconnector.addEntry(amino, pdbid);	
				}
				return true;		
			}catch(SQLException e){
				e.printStackTrace();
				return false;
			}
		} else {
			System.out.println("PDB already exists in DB");
			return false;
		}
	}
		
		public boolean addDSSPEntry(DSSPEntry entry){			
			int pdbDBID = pdbExist(entry.getID()+entry.getChainID()+entry.getChainIDNum());
			int aminoid;
			ResultSet res;
			
			if(pdbDBID != -1 && dsspExist(entry.getID()+entry.getChainID()+entry.getChainIDNum()) == -1){
				try{
					PreparedStatement stmt = connection.createStatement(setDSSPEntry);
					PreparedStatement stmt2 = connection.createStatement("select id from aminoacid where res_index = ? and pdb_id = ?");
					
					for (int i = 0; i < entry.getLength(); i++) {
						
						stmt2.setInt(1, entry.getResIndex()[i]);
						stmt2.setInt(2, pdbDBID);
						res = stmt2.executeQuery();
						res.first();
						aminoid = res.getInt(aafields[0]);
						
						stmt.setString(1, entry.getSecondaryStructure()[i].toString());
						stmt.setInt(2, entry.getAccesability()[i]);
						stmt.setDouble(3, entry.getPhi()[i]);
						stmt.setDouble(4, entry.getPsi()[i]);
						stmt.setDouble(5, entry.getCaTrace()[i][0]);
						stmt.setDouble(6, entry.getCaTrace()[i][1]);
						stmt.setDouble(7, entry.getCaTrace()[i][2]);
						stmt.setInt(8, aminoid);
						stmt.setInt(9, pdbDBID);
						stmt.execute();
					}
					return true;
					
				}catch(SQLException e){
					e.printStackTrace();
					return false;
				}
			} else {
				System.out.println("Corresponding PDB entry doesn't exist or DSSP already imported.");
				return false;
			}
		}
}
