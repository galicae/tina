package db.mysql;

import java.sql.PreparedStatement;
import java.sql.ResultSet;
import java.sql.SQLException;
import java.sql.Statement;
import java.util.ArrayList;
import java.util.List;

import bioinfo.proteins.AminoAcid;
import bioinfo.proteins.Atom;
import bioinfo.proteins.PDBEntry;

public class PDBConnector extends MysqlWrapper{

	private static final String tablename = "pdb";
	
	private static final String[] pdbfields = {"id","pdb_id","chain","length"};
	private static final String[] aafields = {"id","name","res_index","numberofAtom","pdb_id"};
	private static final String[] atomfields = {"id","type","x","y","z","aminoacid_id"};
	private static final String setEntry = "insert into "+tablename+" ("+pdbfields[1]+","+pdbfields[2]+","+pdbfields[3]+") values (?,?,?)";	
	private static final String getById = "select * from pdb join aminoacid on pdb.id = aminoacid.pdb_id join atom on aminoacid.id = atom.aminoacid_id" +
			" where pdb.pdb_id = ?";
	
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
		PreparedStatement stmt = connection.createStatement(getById);
		List<Atom> atoms = new ArrayList<Atom>();
		List<AminoAcid> aminos = new ArrayList<AminoAcid>();
		
		try{
			stmt.setString(1, id);
			ResultSet res = stmt.executeQuery();
			
			int atomquantity;
			double[] pos_temp = new double[3];
					
			while(res.next()){
				//read out aminos and corresponding atoms
				atomquantity = res.getInt(aafields[3]);
				pos_temp[0] = res.getDouble(atomfields[2]);
				pos_temp[1] = res.getDouble(atomfields[3]);
				pos_temp[2] = res.getDouble(atomfields[4]);
				atoms.add(new Atom(res.getString(atomfields[1]),pos_temp));
				
				for (int j = 0; j < atomquantity-1; j++) {
					pos_temp[0] = res.getDouble(atomfields[2]);
					pos_temp[1] = res.getDouble(atomfields[3]);
					pos_temp[2] = res.getDouble(atomfields[4]);
					atoms.add(new Atom(res.getString(atomfields[1]),pos_temp));
					res.next();
				}
				aminos.add(new AminoAcid(aafields[1],res.getInt(aafields[2]),atoms.toArray(new Atom[atomquantity])));
				atoms.clear();
			}
			return new PDBEntry(id,aminos.toArray(new AminoAcid[aminos.size()]));
		}catch(SQLException e){
			e.printStackTrace();
			return null;
		}
	}
	
	public boolean pdbExist(String id){
		Statement stmt = connection.createStatement();
		try{
			ResultSet res = stmt.executeQuery("Select id from "+tablename+" where pdb_id = "+id);
			if(res.first()){
				return true;
			}else{
				return false;
			}
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
			stmt.setInt(3, entry.length());
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
