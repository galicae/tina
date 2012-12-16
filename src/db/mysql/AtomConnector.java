package db.mysql;

import java.sql.PreparedStatement;
import java.sql.ResultSet;
import java.sql.SQLException;
import java.sql.Statement;
import java.util.ArrayList;
import java.util.List;

import bioinfo.proteins.PDBEntry;

public class AtomConnector extends MysqlWrapper{

	private static final String tablename = "atom";
	private static final String[] fields = {"id","type","x","y","z","aminoacid_id"};
	private static final String getById = "select data from atom where id = ?";
	private static final String setEntry = "insert into atom values (?,?)";
	
	public AtomConnector(MysqlDBConnection connection) {
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
			stmt.setString(1,entry.getID());
			stmt.setObject(2, entry);
			
			return stmt.execute();
		}catch(SQLException e){
			e.printStackTrace();
			return false;
		}
	}
	
	
	
	

}
