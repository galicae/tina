package stat.evaltree;

import java.util.Iterator;
import java.util.List;

public abstract class AEvalSupervisor<T> implements IEvalSupervisor<T>{
	
	protected IEvalVertix<T> tree;
	private int totalNumber;
	
	public AEvalSupervisor(){
		generateTree();
	}
	
	public int getTotalNumber(){
		return totalNumber;
	}
	
	public void fillTree(List<T> data){
		Iterator<T> iterator = data.iterator();
		while(iterator.hasNext()){
			tree.evalData(iterator.next());
		}
		totalNumber = data.size();
	}
	
	public void printTree(){
		System.out.println(tree.toString());
	}
	
	public IEvalVertix<T> getByPath(List<String> path){
		if(path.size() == 1){
			if(tree.getLabel().equals(path.get(0))){
				return tree;
			}else{
				return null;
			}
		}else{
			path.remove(0);
			return tree.getByPath(path);
		}
	}
	
	public List<IEvalVertix<T>> getByLabel(String label){
		return tree.getByLabel(label);
	}
	
	public void insertData(T data){
		tree.evalData(data);
	}

}
