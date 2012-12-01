package bioinfo.proteins;

import java.io.BufferedReader;
import java.io.File;
import java.util.ArrayList;
import java.util.List;

/**
 * PDBFileReader read PDBFile and returns it as internal PDBEntry
 * @author andreseitz
 *
 */
public class PDBFileReader {

	private String folder = null;
	private List<String> files = null;;
	
	/**
	 * standard constructor
	 */
	public PDBFileReader(){
		
	}
	
	/**
	 * constructor for forlder-parsing
	 * @param folder
	 */
	public PDBFileReader(String folder){
		this.folder = checkfolder(folder);
	}
	
	public PDBEntry readPDBFromFile(String filename){
		//TODO
		return null;
	}
	
	public List<PDBEntry> readPdbFolder(){
		//TODO
		return null;
	}
	
	public List<PDBEntry> readPdbFolder(String folder){
		//TODO
		return null;
	}
	
	public PDBEntry readFromFolderById(String id){
		//TODO
		return null;
	}
	
	public void initSequentialFolderRead(){
		//TODO
	}
	
	public PDBEntry nextPdb(){
		//TODO
		return  null;
	}
	
	

	public String getFolder() {
		return folder;
	}

	public void setFolder(String folder) {
		this.folder = checkfolder(folder);
	}

	public List<String> getFiles() {
		return files;
	}

	public boolean setFiles(List<String> files) {
		if(folder == null){
			return false;
		}
		File[] existing = new File(folder).listFiles();
		boolean found = false;
		for(String file: files){
			found = false;
			for(int i = 0; i != existing.length; i++){
				if(file.equals(existing[i].getName())){
					found = true;
					break;
				}
			}
			if(!found){
				return false;
			}
		}
		this.files = files;
		return true;
	}
	
	private String checkfolder(String folder){
		if(folder.endsWith("/")){
			return folder;
		}else{
			return folder+"/";
		}
	}
	
	private PDBEntry parseEntry(BufferedReader br, String pdbId){
		List<AminoAcid> aminoacids = new ArrayList<AminoAcid>();
		try{
			char chain = pdbId.charAt(4);
			String line;
			List<Atom> atoms = new ArrayList<Atom>();
			
			String name;
			String resName;
			char chainId;
			int resSeq;
			double[] coord = new double[3];
			
			int lastResSeq = -1;
			String lastResName = "";
			
			
			
			while((line = br.readLine()) != null){
				if(line.startsWith("ATOM")){
					chainId = line.charAt(21);
					if(chainId != chain){
						continue;
					}
					resSeq = Integer.parseInt(line.substring(22,26).trim());
					name = line.substring(12,16).trim();
					resName = line.substring(17,20).trim();
					coord[0] = Double.parseDouble(line.substring(30,38).trim());
					coord[1] = Double.parseDouble(line.substring(38,46).trim());
					coord[2] = Double.parseDouble(line.substring(46,54).trim());
					
					if(lastResSeq != resSeq){
						aminoacids.add(new AminoAcid(AminoAcidName.getAAFromTLC(resName),atoms.toArray(new Atom[atoms.size()])));
						atoms.clear();
					}
					lastResSeq = resSeq;
					lastResName = resName;
					atoms.add(new Atom(name,coord));
				}
				aminoacids.add(new AminoAcid(AminoAcidName.getAAFromTLC(lastResName),atoms.toArray(new Atom[atoms.size()])));
			}
		}catch(Exception e){
			e.printStackTrace();
		}
		return new PDBEntry(pdbId, aminoacids);
	}
	
	

	
}
