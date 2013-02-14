package stat.evaltree;

import java.util.ArrayList;
import java.util.List;

public abstract class AEvalLeaf<T> extends AEvalVertix<T>{

	private List<T> data = new ArrayList<T>();
	
	public AEvalLeaf(String label, String description) {
		super(label, description);
	}
	
	public int getNumber(){
		return data.size();
	}
	
	public boolean isLeaf(){
		return true;
	}
	
	public List<IEvalVertix<T>> getChildren(){
		return null;
	}
	
	@Override
	public List<IEvalVertix<T>> getByLabel(String label){
		List<IEvalVertix<T>> result = new ArrayList<IEvalVertix<T>>();
		if(this.label.equals(label)){
			result.add(this);
		}
		return result;
	}
	
	public IEvalVertix<T> getByPath(List<String> path){
		return null;
	}
	
	public String toString(){
		return "\n"+label+"\t"+getNumber();
	}

	public List<T> getData() {
		return data;
	}

	public void setData(List<T> data) {
		this.data = data;
	}
	
	public void addData(T data){
		this.data.add(data);
	}
}
