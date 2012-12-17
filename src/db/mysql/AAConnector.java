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

public class AAConnector extends MysqlWrapper{

	private static final String tablename = "aminoacid";
	private static final String[] fields = {"id","name","res_index","numberofAtom","pdb_id"};
	private static final String getById = "select data from aminoacid where id = ?";
	private static final String setEntry = "insert into aminoacid ("+fields[1]+","+fields[2]+","+fields[3]+","+fields[4]+") values (?,?,?,?)";
	
	public AAConnector(MysqlDBConnection connection) {
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
			lastid = res.getInt(0);
		} catch (SQLException e) {
			e.printStackTrace();
			return 0;
		}	
		return lastid;
	}
	
	public boolean addEntry(AminoAcid entry, int pdbid){
		AtomConnector atomconnector = new AtomConnector(this.connection);
		PreparedStatement stmt = connection.createStatement(setEntry);
				
		try{
			Atom atom;
			
			stmt.setString(1, entry.toString());
			stmt.setInt(2, entry.getResIndex());
			stmt.setInt(3, entry.getAtomNumber());
			stmt.setInt(4, pdbid);
			int aminoid = getLastId();
					
			for (int i = 0; i < entry.getAtomNumber(); i++) {
				//insert atoms
				atom = entry.getAtom(i);
				atomconnector.addEntry(atom, aminoid);	
			}
			return stmt.execute();
		}catch(SQLException e){
			e.printStackTrace();
			return false;
		}
	}
	
	
	
	

}
