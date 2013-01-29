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

public class PDBConnector extends MysqlWrapper{

	private static final String tablename = "pdb";
	
	private static final String[] pdbfields = {"id","pdb_id","chainID", "chainIDNum","length"};
	private static final String[] dsspfields = {"id","ss","sa", "phi","psi","x_ca","y_ca","z_ca","seq_index","aminoacid_id","pdb_id"};
	private static final String[] aafields = {"id","name","res_index","numberofAtom","pdb_id"};
	private static final String[] atomfields = {"id","type","x","y","z","aminoacid_id"};
	private static final String setEntry = "insert into "+tablename+" ("+pdbfields[1]+","+pdbfields[2]+","+pdbfields[3]+","+pdbfields[4]+") values (?,?,?,?)";	
	private static final String getPDBById = "select * from pdb join aminoacid on pdb.id = aminoacid.pdb_id join atom on aminoacid.id = atom.aminoacid_id" +
			" where pdb.pdb_id = ? and pdb.chainID = ? and pdb.chainIDNum = ?";
	private static final String getDSSPById = "select dssp.*,aminoacid.name from dssp,aminoacid join pdb on pdb.id = dssp.pdb_id join aminoacid.id = dssp.aminoacid_id" +
			" where pdb.pdb_id = ? and pdb.chainID = ? and pdb.chainIDNum = ?";
	
	public PDBConnector(MysqlDBConnection connection) {
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
		int chainIDNum=Integer.valueOf(id.substring(5, 7));
		String pdbid = id.substring(0, 4);
		
		try{
			stmt.setString(1, pdbid);
			stmt.setString(2, chainID+"");
			stmt.setInt(3, chainIDNum);
			ResultSet res = stmt.executeQuery();
			
			int atomquantity;
			double[] pos_temp;
			
			if(pdbExist(id)){
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
		List<SecStructEight> secStruct = new ArrayList<SecStructEight>();
		List<double[]> caTrace = new ArrayList<double[]>();
		List<Integer> sa = new ArrayList<Integer>();
		List<Double> phi = new ArrayList<Double>();
		List<Double> psi = new ArrayList<Double>();
		
		char chainID = id.charAt(5);
		int chainIDNum=Integer.valueOf(id.substring(5, 7));
		id = id.substring(0, 4);

		try{
			stmt.setString(1, id);
			stmt.setString(2, chainID+"");
			stmt.setInt(3, chainIDNum);
			ResultSet res = stmt.executeQuery();

			if(pdbExist(id)){		
				while(res.next()){
					//read out dssp
					aminos.add(AminoAcidName.getAAFromOLC(res.getString(aafields[1])));
					secStruct.add(SecStructEight.getSSFromChar(res.getString(dsspfields[1]).charAt(0)));
					caTrace.add(new double[] {res.getDouble(dsspfields[5]),res.getDouble(dsspfields[6]),res.getDouble(dsspfields[7])});
					sa.add(res.getInt(dsspfields[2]));
					phi.add(res.getDouble(dsspfields[3]));
					phi.add(res.getDouble(dsspfields[4]));
					
				}
				return new DSSPEntry(id,aminos,secStruct,sa,phi,psi,caTrace);
			} else {
				return null;
			}
		}catch(SQLException e){
			e.printStackTrace();
			return null;
		}
	}
	
	public boolean pdbExist(String id){
		String query = "select * from "+tablename+" where "+pdbfields[1]+" = ? and "+pdbfields[2]+" = ? and "+pdbfields[3]+" = ?";
		PreparedStatement stmt = connection.createStatement(query);
		
		String pdbid = id.substring(0,4);
		String chain = id.substring(4,5);
		int chainIdNum = Integer.parseInt(id.substring(5,7));

		try{
			stmt.setString(1, pdbid);
			stmt.setString(2, chain);
			stmt.setInt(3, chainIdNum);
			
			ResultSet res = stmt.executeQuery();
			
			if(res.first()){
				return true;
			}
			return false;
		}catch(SQLException e){
			e.printStackTrace();
			return false;
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
	
	public boolean addEntry(PDBEntry entry){
		AAConnector aaconnector = new AAConnector(this.connection);
		PreparedStatement stmt = connection.createStatement(setEntry);
		
		try{
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
	}
	
	
	
	

}
