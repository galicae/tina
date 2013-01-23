package opt;

import java.util.List;

public class Setting {
	
	private String optionIdentifier;
	private boolean required;
	private Object standard;
	private ParameterType type;
	private boolean validate;
	private List<ValidationRule> validations;
	private String outputIdentifier;
	
	public String getOptionIdentifier() {
		return optionIdentifier;
	}
	public void setOptionIdentifier(String optionIdentifier) {
		this.optionIdentifier = optionIdentifier;
	}
	public boolean isRequired() {
		return required;
	}
	public void setRequired(boolean required) {
		this.required = required;
	}
	public ParameterType getType() {
		return type;
	}
	public void setType(ParameterType type) {
		this.type = type;
	}
	public boolean isValidate() {
		return validate;
	}
	public void setValidate(boolean validate) {
		this.validate = validate;
	}
	public List<ValidationRule> getValidations() {
		return validations;
	}
	public void setValidations(List<ValidationRule> validations) {
		this.validations = validations;
	}
	public String getOutputIdentifier() {
		return outputIdentifier;
	}
	public void setOutputIdentifier(String outputIdentifier) {
		this.outputIdentifier = outputIdentifier;
	}
	public Object getStandard() {
		return standard;
	}
	public void setStandard(Object standard) {
		this.standard = standard;
	}
	public Setting(String optionIdentifier, String outputIdentifier,boolean required, Object standard,
			ParameterType type, boolean validate, List<ValidationRule> validations) {
		this.optionIdentifier = optionIdentifier;
		this.required = required;
		this.type = type;
		this.validate = validate;
		this.validations = validations;
		this.outputIdentifier = outputIdentifier;
		this.standard = standard;
	}
	
	

}
