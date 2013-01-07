package db.mysql;

import java.sql.PreparedStatement;
import java.sql.ResultSet;
import java.sql.SQLException;
import java.sql.Statement;
import java.util.ArrayList;
import java.util.List;

import bioinfo.proteins.Atom;
import bioinfo.proteins.PDBEntry;

public class AtomConnector extends MysqlWrapper{

	private static final String tablename = "atom";
	private static final String[] fields = {"id","type","x","y","z","aminoacid_id"};
	private static final String getById = "select data from atom where id = ?";
	private static final String setEntry = "insert into atom ("+fields[1]+","+fields[2]+","+fields[3]+","+fields[4]+","+fields[5]+") values (?,?,?,?,?)";
	
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
	
	public boolean addEntry(Atom entry, int aminoid){
		PreparedStatement stmt = connection.createStatement(setEntry);
				
		try{
			double[] pos = entry.getPosition();
			stmt.setString(1, entry.getType().toString());
			stmt.setDouble(2, pos[0]);
			stmt.setDouble(3, pos[1]);
			stmt.setDouble(4, pos[2]);
			stmt.setInt(5, aminoid);
			return stmt.execute();
		}catch(SQLException e){
			e.printStackTrace();
			return false;
		}
	}
}
