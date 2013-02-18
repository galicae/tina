package bioinfo.energy.potential;

import java.io.Writer;
import java.util.HashMap;

public class PotentialDimension<A extends Number> {
	
	private PotentialDimension<A>[] subDimensions = null;
	private A value;
	private final int level;
	
	public PotentialDimension(int[] size, A initVal){
		this(size, 0, initVal);
	}
	
	@SuppressWarnings("unchecked")
	private PotentialDimension(int[] size, int level, A initVal){
		this.level = level;
		if(size.length >= level+1){
			//has subDimensions
			this.subDimensions = (PotentialDimension<A>[])new PotentialDimension[size[level]];
			for(int i = 0; i != this.subDimensions.length; i++){
				this.subDimensions[i] = new PotentialDimension<A>(size, level+1, initVal);
			}
		}else{
			// has no subDimensions
			this.subDimensions = null;
			this.value = initVal;
		}
	}
	
	public boolean hasSubDimensions(){
		if(subDimensions != null && subDimensions.length > 0){
			return true;
		}
		return false;
	}
	
	public PotentialDimension<A> getByAddress(int[] address){
		if(address.length == level){
			return this;
		}else{
			if(this.hasSubDimensions()){
				if(this.subDimensions.length >address[level]){
					return subDimensions[address[level]].getByAddress(address);
				}else{
					return null;
				}
			}else{
				return null;
			}
		}
	}
		
	public boolean setValue(int[] address, A value){
		if(address.length == level){
			if(this.hasSubDimensions()){
				return false;
			}else{
				this.value = value;
				return true;
			}
		}else{
			if(this.subDimensions.length >address[level]){
				return this.subDimensions[address[level]].setValue(address, value);
			}else{
				return false;
			}
		}
	}
	
	public boolean writeAsXML(Writer stream){
		try{
			if(this.hasSubDimensions()){
				stream.append("<dimension level=\""+this.level+"\">");
				stream.append("\n");
				for(int i = 0; i != this.subDimensions.length; i++){
					this.subDimensions[i].writeAsXML(stream);
				}
				stream.append("</dimension>");
				stream.append("\n");
			}else{
				if(value != null){
					stream.append("<data value=\""+value.toString()+"\"/>");
				}else{
					stream.append("<data value=\"0\"/>");
				}

				stream.append("\n");
			}
			return true;
		}catch(Exception e){
			e.printStackTrace();
			return false;
		}
	}

	public HashMap<Integer,Integer> getSize(){
		HashMap<Integer,Integer> size = new HashMap<Integer,Integer>();
		if(hasSubDimensions()){
			for(int i = 0; i != this.subDimensions.length; i++){
				size.putAll(this.subDimensions[i].getSize());
			}
			size.put(this.level, this.subDimensions.length);
		}
		return size;
	}
	
	public A getValue(){
		if(!hasSubDimensions()){
			return value;
		}else{
			System.err.println("out of bounds");
			return null;
		}
	}
	

}
