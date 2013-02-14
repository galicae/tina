package stat.evaltree;

import java.util.List;

public interface IEvalSupervisor<T> {

	public void generateTree();
	public void fillTree(List<T> data);
	public void insertData(T data);
	public void printTree();
	public IEvalVertix<T> getByPath(List<String> path);
	public List<IEvalVertix<T>> getByLabel(String label);

	
}
