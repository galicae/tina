package bioinfo.proteins;

import bioinfo.alignment.Alignable;

public class SSCCEntry implements Alignable {

	private SSCCLine[] lines;
	private int length;
	private final String id;
	
	public SSCCEntry(String id,SSCCLine[] lines){
		this.id = id;
		this.lines = lines;
		this.length = lines.length;
	}

	@Override
	public Object getComp(int i) {
		return this.lines[i];
	}

	@Override
	public int length() {
		return this.length;
	}

	@Override
	public String getId() {
		return this.id;
	}

}
