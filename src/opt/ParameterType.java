package opt;

public enum ParameterType {
	Integer(new Integer(0).getClass()),
	Double(new Double(0.0).getClass()),
	String(new String("").getClass()),
	Character(new Character((char)0).getClass()),
	Boolean(new Boolean(true).getClass());
	
	private final Object objecttype;
	
	private ParameterType(Object object){
		this.objecttype = object.getClass();
	}
	
	public String classString(){
		for(ParameterType x: ParameterType.values()){
			if(x.equals(this)){
				return this.objecttype.toString();
			}
		}
		return null;
	}
	
	public Object cast(String x) throws Exception{
		switch(this){
			case Boolean:{
				return true;
			}
			case Integer:{
				return java.lang.Integer.parseInt(x);
			}
			case String:{	
				return x;
			}
			case Character:{
				char[] tmp = x.toCharArray();
				if(tmp.length > 1){
					throw new Exception("\""+x+"\""+" is no character");
				}
				return tmp[0];
			}
			case Double:{
				return java.lang.Double.parseDouble(x);
			}
			default:{
				return null;
			}
			
		}
	}
	
	

}
