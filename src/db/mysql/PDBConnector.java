package db.mysql;

import java.sql.PreparedStatement;
import java.sql.ResultSet;
import java.sql.SQLException;
import java.sql.Statement;
import java.util.ArrayList;
import java.util.List;

import bioinfo.proteins.PDBEntry;

public class PDBConnector extends MysqlWrapper{

	private static final String tablename = "pdb";
	private static final String[] fields = {"id","pdb_id","chain","length"};
	private static final String getById = "select data from pdb where id = ?";
	private static final String setEntry = "insert into pdb ("+fields[1]+","+fields[2]+","+fields[3]+") values (?,?,?)";
	
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
	
	public boolean addEntry(PDBEntry entry){
		PreparedStatement stmt = connection.createStatement(setEntry);
		try{
			//insert pdbentry
			stmt.setString(1,entry.getID());
			stmt.setString(2,String.valueOf(entry.getChainID()));
			stmt.setInt(3, entry.length());;	
			return stmt.execute();
		}catch(SQLException e){
			e.printStackTrace();
			return false;
		}
	}
	
	
	
	

}
