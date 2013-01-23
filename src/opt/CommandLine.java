package opt;

import java.util.HashMap;
import java.util.List;

public class CommandLine {
	
	private String[] commands;
	private List<Setting> settings;
	private CommandParser parser;
	private HashMap<String,Object> parsedComs;
	
	
	public String[] getCommands() {
		return commands;
	}
	public void setCommands(String[] commands) {
		this.commands = commands;
	}
	public List<Setting> getSettings() {
		return settings;
	}
	public void setSettings(List<Setting> settings) {
		this.settings = settings;
	}
	public HashMap<String, Object> getParsedComs() {
		return parsedComs;
	}
	public void setParsedComs(HashMap<String, Object> parsedComs) {
		this.parsedComs = parsedComs;
	}
	public CommandParser getParser() {
		return parser;
	}
	public void setParser(CommandParser parser) {
		this.parser = parser;
	}
	
	public CommandLine(List<Setting> settings,String[] commands) {
		this.commands = commands;
		this.settings = settings;
		this.parser = new CommandParser(settings,commands);
		try{
			this.parsedComs = parser.parse();
		} catch(Exception e){
			e.printStackTrace();
			this.parsedComs = null;
		}
	}
	
	
	public CommandLine(List<Setting> settings) {
		this.settings = settings;
	}
	
	public void parse(String[] commands){
		this.parser = new CommandParser(settings,commands);
		try{
			this.parsedComs = parser.parse();
		} catch(Exception e){
			e.printStackTrace();
			this.parsedComs = null;
		}
	}
	
	
	
	
	

}
