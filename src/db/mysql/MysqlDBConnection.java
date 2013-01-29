package db.mysql;

import java.sql.PreparedStatement;
import java.sql.Statement;


public interface MysqlDBConnection {
		
		
	public boolean isAlive();
	
	public Statement createStatement();
	
	public PreparedStatement createStatement(String preparation);
	
	public boolean containsTable(String tablename);
	
	

}
