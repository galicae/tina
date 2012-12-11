package db.mysql;

import java.sql.ResultSet;
import java.sql.SQLException;
import java.sql.Statement;

public abstract class MysqlWrapper {
	
	/**
	 * class variable holding db connection which should be initialized only once to minimize network traffic
	 */
	protected MysqlDBConnection connection = null;

	
	/**
	 * 
	 * @param connection to db, should only be initialized once for all wrappers to minimalize network traffic
	 */
	public MysqlWrapper(MysqlDBConnection connection){
		this.connection = connection;
	}
	
	/**
	 * 
	 * @return size of the table
	 */
	public int size() {
		Statement stmt = connection.createStatement();
		try{
			ResultSet res = stmt.executeQuery("select count(*) as size from "+getTablename());
			return res.getInt("size");
		} catch(SQLException e){
			e.printStackTrace();
			return -1;
		}
	}
	
	/**
	 * 
	 * @return tablename of specific tableWrapper
	 * tablename should be held in a static final class variable
	 */
	abstract String getTablename();
	
	/**
	 * 
	 * @return Array of field names of specific tableWrapper
	 * fieldnames should be held in a static final class variable
	 */
	abstract String[] getFields();
	
	/**
	 * 
	 * @return true if connection is up and running and table exists
	 *
	 */
	public boolean isAccessible(){
		return connection.isAlive() && connection.containsTable(getTablename());
	}
	
}
