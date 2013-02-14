package stat.evaltree;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;

public abstract class AEvalNode<T> extends AEvalVertix<T>{
	
	protected List<IEvalVertix<T>> children;
	protected HashMap<String,Integer> classifications;
	
	public AEvalNode(String label, String description, List<IEvalVertix<T>> children, HashMap<String,Integer> classifications) {
		super(label, description);
		this.children = children;
		this.classifications = classifications;
	}
	
	public int getNumber(){
		int number = 0;
		Iterator<IEvalVertix<T>> children = this.children.iterator();
		while(children.hasNext()){
			number += children.next().getNumber();
		}
		return number;
	}
	
	public boolean isLeaf(){
		return false;
	}
	
	public List<IEvalVertix<T>> getChildren(){
		return children;
	}
	
	public List<IEvalVertix<T>> getByLabel(String label){
		List<IEvalVertix<T>> result = new ArrayList<IEvalVertix<T>>();
		if(this.label.equals(label)){
			result.add(this);
		}
		Iterator<IEvalVertix<T>> children = this.children.iterator();
		Iterator<IEvalVertix<T>> childsChildren;
		while(children.hasNext()){
			childsChildren = children.next().getByLabel(label).iterator();
			while(childsChildren.hasNext()){
				result.add(childsChildren.next());
			}
			
		}
		return result;
	}
	
	public IEvalVertix<T> getByPath(List<String> path){
		if(path.size() == 1){
			Iterator<IEvalVertix<T>> iterator = children.iterator();
			IEvalVertix<T> tmp;
			while(iterator.hasNext()){
				if((tmp = iterator.next()).getLabel().equals(path.get(0))){
					return tmp;
				}
			}
			return null;
		}else{
			Iterator<IEvalVertix<T>> iterator = children.iterator();
			IEvalVertix<T> tmp;
			IEvalVertix<T> result;
			while(iterator.hasNext()){
				if((tmp = iterator.next()).getLabel().equals(path.get(0))){
					path.remove(0);
					if((result = tmp.getByPath(path)) != null){
						return result;
					}
				}
			}
			return null;
		}
	}
	
	public String toString(){
		String result = "\n"+label+"\t"+getNumber();
		Iterator<IEvalVertix<T>> iterator = children.iterator();
		while(iterator.hasNext()){
			result += iterator.next().toString().replaceAll("\n", "\n\t");
		}
		return result;
	}

	public void setChildren(List<IEvalVertix<T>> children) {
		this.children = children;
	}
	
	

}
