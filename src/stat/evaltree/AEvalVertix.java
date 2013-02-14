package stat.evaltree;

public abstract class AEvalVertix<T> implements IEvalVertix<T>{

	String label;
	String description;
	
	public AEvalVertix(String label, String description){
		this.label = label;
		this.description = description;
	}
	
	public String getDescription(){
		return description;
	}
	
	public String getLabel(){
		return label;
	}
	
	
}
