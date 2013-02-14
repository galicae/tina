package stat.evaltree;

import java.util.List;

public interface IEvalVertix<T> {
	
	public String getLabel();
	public String getDescription();
	public void evalData(T data);
	public int getNumber();
	public boolean isLeaf();
	public List<IEvalVertix<T>> getChildren();
	public List<IEvalVertix<T>> getByLabel(String label);
	public IEvalVertix<T> getByPath(List<String> path);
	public String toString();
	

}
