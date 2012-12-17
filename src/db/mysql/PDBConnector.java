package db.mysql;

import java.sql.PreparedStatement;
import java.sql.ResultSet;
import java.sql.SQLException;
import java.sql.Statement;
import java.util.ArrayList;
import java.util.List;

import bioinfo.proteins.AminoAcid;
import bioinfo.proteins.PDBEntry;

public class PDBConnector extends MysqlWrapper{

	private static final String tablename = "pdb";
	
	private static final String[] fields = {"id","pdb_id","chain","length"};
	private static final String setEntry = "insert into "+tablename+" ("+fields[1]+","+fields[2]+","+fields[3]+") values (?,?,?)";	
	private static final String getById = "select * from pdb where pdb_id = ?";
	
	public PDBConnector(MysqlDBConnection connection) {
		super(connection);
	}
	
	@Override
	String getTablename() {
		return tablename;
	}

	@Override
	String[] getFields() {
		return fields;
	}
	
	public PDBEntry getPDB(String id){
		PreparedStatement stmt = connection.createStatement(getById);
		try{
			stmt.setString(1, id);
			ResultSet res = stmt.executeQuery();
			if(res.first()){
				return (PDBEntry)res.getObject(fields[1]);
			}else{
				return null;
			}
		}catch(SQLException e){
			e.printStackTrace();
			return null;
		}
	}
	
	public String[] getIDs(){
		Statement stmt = connection.createStatement();
		try {
			ResultSet res = stmt.executeQuery("select id from "+getTablename());
			List<String> ids = new ArrayList<String>();
			while(res.next()){
				ids.add(res.getString(fields[0]));
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
