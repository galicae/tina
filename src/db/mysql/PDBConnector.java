package db.mysql;

import java.sql.ResultSet;
import java.sql.SQLException;
import java.sql.Statement;

import com.mysql.jdbc.Connection;

public class PDBConnector extends MysqlWrapper{

	private static final String tablename = "pdb";
	private static final String[] fields = {"id","data"};
	
	public PDBConnector(MysqlDBConnection connection) {
		super(connection);
	}

	@Override
	int size() {
		Statement stmt = connection.createStatement();
		try{
			ResultSet res = stmt.executeQuery("select count(*) as size from "+tablename);
			return res.getInt("size");
		} catch(SQLException e){
			e.printStackTrace();
			return -1;
		}
	}

	@Override
	String getTablename() {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	String[] getFields() {
		// TODO Auto-generated method stub
		return null;
	}
	
	

}
