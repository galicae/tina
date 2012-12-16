package bioinfo.proteins;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.List;

/**
 * PDBFileReader read PDBFile and returns it as internal PDBEntry
 * @author andreseitz
 *
 */
public class PDBFileReader {

	private String folder = null;
	private List<String> files = null;
	private int sequentialReadIndex = 0;
	
	/**
	 * standard constructor creating "empty" PDBFileReader
	 */
	public PDBFileReader(){
		
	}
	
	/**
	 * checked folder will end with "/"
	 * @param folder which ends with "/"
	 */
	public PDBFileReader(String folder){
		this.folder = checkfolder(folder);
	}
	
	/**
	 * 
	 * @param filename representing pdb file to read in format: relpath/aaaaA00.pdb
	 * @return
	 */
	public PDBEntry readPDBFromFile(String filename){
		BufferedReader br = null;
		try{
			 br = new BufferedReader(new InputStreamReader(new FileInputStream(filename)));
		}catch(Exception e){
			e.printStackTrace();
		}
		String id = filename.split("/")[filename.split("/").length-1].split("\\.")[0];
		return parseEntry(br, id);
	}
	
	/**
	 * @return List of PDBEntries in folder, null if no pdbs are in folder or folder is not set
	 */
	public List<PDBEntry> readPdbFolder(){
		List<PDBEntry> result = new ArrayList<PDBEntry>();
		BufferedReader br = null;
		String id;
		try{
			if(folder == null){
				return null;
			}
			File[] files = new File(folder).listFiles();
			if(files.length == 0){
				return null;
			}
			for(File file : files){
				br = new BufferedReader(new InputStreamReader(new FileInputStream(folder+file)));
				id = file.getName().split("\\.")[0];
				result.add(parseEntry(br,id));
			}
			return result;
		}catch(Exception e){
			e.printStackTrace();
		}
		return null;
	}
	
	/**
	 * 
	 * @param folder to read pdbs from, folder will be set as folder-class-variable
	 * @return List of PDBEntries in folder, null if no pdbs are in folder or folder is not set
	 */
	public List<PDBEntry> readPdbFolder(String folder){
		this.folder = checkfolder(folder);
		List<PDBEntry> result = new ArrayList<PDBEntry>();
		BufferedReader br = null;
		String id;
		try{
			if(folder == null){
				return null;
			}
			File[] files = new File(folder).listFiles();
			if(files.length == 0){
				return null;
			}
			for(File file : files){
				br = new BufferedReader(new InputStreamReader(new FileInputStream(folder+file)));
				id = file.getName().split("\\.")[0];
				result.add(parseEntry(br,id));
			}
			return result;
		}catch(Exception e){
			e.printStackTrace();
		}
		return null;
	}
	
	/**
	 * 
	 * @param id of pdb to read from
	 * @return PDBentry representing the pdb of given id, null if folder is not set or certain pdb does not exist
	 */
	public PDBEntry readFromFolderById(String id){
		if(folder == null){
			return null;
		}
		// huberste: PDBFiles normally are named without the ChainID!
		String filename = folder+id.substring(0, 4).toUpperCase()+".pdb"; 
		BufferedReader br = null;
		try{
			br = new BufferedReader(new InputStreamReader(new FileInputStream(filename)));
		}catch(Exception e){
			return null;
		}
		return parseEntry(br,id);
	}
	
	/**
	 * initialise PDBEntry list to read it out sequentially
	 */
	public void initSequentialFolderRead(){
		if(folder == null){
			return;
		}
		File[] files = new File(folder).listFiles();
		List<String> names = new ArrayList<String>();
		for(File file: files){
			names.add(file.getName());
		}
		this.files = names;
		sequentialReadIndex = 0;
	}
	
	/**
	 * 
	 * @return next PDBEntry from initialisedList
	 */
	public PDBEntry nextPdb(){
		if(files.size() > sequentialReadIndex){
			String file = files.get(sequentialReadIndex);
			BufferedReader br;
			try{
				br = new BufferedReader(new InputStreamReader(new FileInputStream(folder+file)));
				String id = file.split("\\.")[0];
				sequentialReadIndex++;
				return parseEntry(br, id);
			}catch(Exception e){
				e.printStackTrace();
			}
			return null;
		}else{
			return null;
		}
		
	}
	
	/**
	 * 
	 * @return folder to read from
	 */
	public String getFolder() {
		return folder;
	}

	/**
	 * 
	 * @param folder to read from
	 */
	public void setFolder(String folder) {
		this.folder = checkfolder(folder);
	}

	/**
	 * 
	 * @return list of files being set for read out 
	 */
	public List<String> getFiles() {
		return files;
	}

	/**
	 * 
	 * @param files which shall be read from folder if folder != null or from
	 * anywhere in the file system if folder == null
	 * @return true if files were set correctly or false if setFiles failed
	 */
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
	
	/**
	 * 
	 * @param folder which should be checked for validity
	 * @return
	 */
	private String checkfolder(String folder){
		if(folder.endsWith("/")){
			return folder;
		}else{
			return folder+"/";
		}
	}
	
	/**
	 * 
	 * @param br BufferedReader to read from
	 * @param pdbId pdbId of format aaaaA00
	 * @return PDbEntry representing the read file
	 */
	private PDBEntry parseEntry(BufferedReader br, String pdbId){
		List<AminoAcid> aminoacids = new ArrayList<AminoAcid>();
		try{
			char chain = pdbId.charAt(4);
			String line;
			List<Atom> atoms = new ArrayList<Atom>();
			
			String name;
			String resName = "";
			char chainId;
			int resSeq = 0;
			double[] coord;
			
			int lastResSeq = 0;
			
			while((line = br.readLine()) != null){
				if(line.startsWith("ATOM")){
					chainId = line.charAt(21);
					if(chainId != chain){
						continue;
					}
					resSeq = Integer.parseInt(line.substring(22,26).trim());
					name = line.substring(12,16).trim();
					resName = line.substring(17,20).trim();
					coord = new double[3];
					coord[0] = Double.parseDouble(line.substring(30,38).trim());
					coord[1] = Double.parseDouble(line.substring(38,46).trim());
					coord[2] = Double.parseDouble(line.substring(46,54).trim());
					
					if(lastResSeq != resSeq){
						aminoacids.add(new AminoAcid(AminoAcidName.getAAFromTLC(resName),resSeq,atoms.toArray(new Atom[atoms.size()])));
						atoms.clear();
						lastResSeq = resSeq;
					}
					atoms.add(new Atom(name,coord));
				}
			}
			if(atoms != null && atoms.size() != 0){
				aminoacids.add(new AminoAcid(AminoAcidName.getAAFromTLC(resName),resSeq,atoms.toArray(new Atom[atoms.size()])));
			}
			br.close();
		}catch(Exception e){
			e.printStackTrace();
		}
		return new PDBEntry(pdbId, aminoacids);
	}

}
