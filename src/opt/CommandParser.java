package opt;

import java.util.HashMap;
import java.util.List;

import javax.management.ObjectInstance;

public class CommandParser{
	
	private List<Setting> settings;
	private String[] commandLine;
	
	public CommandParser(List<Setting> settings, String[] commandLine) throws Exception{
		this.settings = settings;
		this.commandLine = commandLine;
	}
	
	public HashMap<String,Object> parse() throws Exception{
		HashMap<String,Object> result = new HashMap<String, Object>();
		for(int pos = 0; pos < commandLine.length; pos++){
			boolean doneFlag = false;
			for(Setting set: settings){
				if(set.getOptionIdentifier().equals(commandLine[pos])){
					switch(set.getType()){
						case Boolean:{
							result.put(set.getOutputIdentifier(), true);
							break;
						}
						default:{
							String value = commandLine[++pos];
							if(set.isValidate()){
								for(ValidationRule validation: set.getValidations()){
									if(!validation.validate(value)){
										throw new Exception("Validation of "+set.getOptionIdentifier()+" failed cause of value "+value);	
									}
								}
								result.put(set.getOutputIdentifier(),set.getType().cast(value));
							}else {
								result.put(set.getOutputIdentifier(),set.getType().cast(value));
							}
						}
					}
					doneFlag = true;
					break;
				}
			}
			if(!doneFlag){
				throw new Exception("Option "+commandLine[pos]+" couldn't be interpreted by parser");

			}
		}
		for(Setting set: settings){
			if(set.isRequired() && !result.containsKey(set.getOutputIdentifier())){
				if(set.getStandard() != null){
					result.put(set.getOutputIdentifier(), set.getStandard());
				} else {
					throw new Exception(set.getOptionIdentifier()+" wasnt determined in commadLine");
				}
			}
		}
		return result;
	}
	
	

}
