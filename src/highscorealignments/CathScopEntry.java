package highscorealignments;

public class CathScopEntry {
	private String id;
	
	//cath info
	private int cath_family;
	private int cath_superfamily;
	private int cath_fold;
	private char cath_clazz;
	
	//scop info
	private int scop_family;
	private int scop_superfamily;
	private int scop_fold;
	private int scop_clazz;	
	
	public CathScopEntry(String id, char cath_c,int cath_fam, int cath_supfam, int cath_fold, int scop_c, int scop_fam, int scop_supfam, int scop_fold){
		this.id = id;
		
		this.cath_family = cath_fam;
		this.cath_superfamily = cath_supfam;
		this.cath_fold = cath_fold;
		this.cath_clazz = cath_c;
		
		this.scop_family = scop_fam;
		this.scop_superfamily = scop_supfam;
		this.scop_fold = scop_fold;
		this.scop_clazz = scop_c;
	}
	public String getID(){
		return this.id;
	}
	
	//cath
	public char getCathClazz(){
		return this.cath_clazz;
	}
	public int getCathFam(){
		return this.cath_family;
	}
	public int getCathSupFam(){
		return this.cath_superfamily;
	}

	public int getCathFold(){
		return this.cath_fold;
	}
	
	//scop
	public int getScopClazz(){
		return this.scop_clazz;
	}
	public int getScopFam(){
		return this.scop_family;
	}
	public int getScopSupFam(){
		return this.scop_superfamily;
	}

	public int getScopFold(){
		return this.scop_fold;
	}
}
