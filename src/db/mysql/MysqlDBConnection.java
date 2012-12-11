package db.mysql;

import java.sql.Connection;
import java.sql.DatabaseMetaData;
import java.sql.DriverManager;
import java.sql.PreparedStatement;
import java.sql.ResultSet;
import java.sql.SQLException;
import java.sql.Statement;
import java.util.Properties;

public class MysqlDBConnection {
	
	private static final String HOST = "mysql2-ext.bio.ifi.lmu.de";
	private static final String USER = "gobig4";
	private static final String PASS = "sHpKnPaS";
	private static final String DB = "gobig4";
	
	private Connection connection = null;
	
	public MysqlDBConnection(){
		Properties connectionProps = new Properties();
	    connectionProps.put("user", USER);
	    connectionProps.put("password", PASS);
	    
	    try {
			connection = DriverManager.getConnection("jdbc:mysql://" + HOST + "/", connectionProps);
		} catch (SQLException e) {
			System.err.println("Mysql connection fail! user:"+USER+" pass:"+PASS+" host:"+HOST+" db:"+DB);
			e.printStackTrace();
		}
	}
	
	public boolean isAlive(){
		try {
			if(connection.isClosed()){
				return false;
			}else{
				return true;
			}
		} catch (SQLException e) {
			System.err.println("Error occured during mysql-dbv live check!");
			return false;
		}
	}
	
	public Statement createStatement(){
		try {
			if(!isAlive()){
				return null;
			}
			return connection.createStatement();
		} catch (SQLException e) {
			e.printStackTrace();
			return null;
		}
	}
	
	public PreparedStatement createStatement(String preparation){
		try {
			if(!isAlive()){
				return null;
			}
			return connection.prepareStatement(preparation);
		} catch (SQLException e) {
			e.printStackTrace();
			return null;
		}
	}
	
	public boolean containsTable(String tablename){
		try{
			DatabaseMetaData meta = connection.getMetaData();
			ResultSet res = meta.getTables(null, null, null, new String[] {"TABLE"});
			while (res.next()){
				if(res.getString("TABLE_NAME").equals(tablename)){
					return true;
				}
			}
			return false;
		} catch(SQLException e){
			e.printStackTrace();
			return false;
		}
	}
	
	

}
