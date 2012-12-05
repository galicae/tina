package bioinfo.alignment;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.List;

import bioinfo.Sequence;
import bioinfo.alignment.gotoh.Gotoh;

/**
 * SequenceAlignmentFileReader read PDBFile and returns it as internal Alignment
 * @author andreseitz
 *
 */
public class SequenceAlignmentFileReader {

	private String filename = null;
	private BufferedReader br;
	
	public SequenceAlignmentFileReader(){
		
	}
	
	public SequenceAlignmentFileReader(String filename){
		this.filename = filename;
	}
	
	public List<SequenceAlignment> readAlignments(){
		List<SequenceAlignment> alignments = new ArrayList<SequenceAlignment>();
		BufferedReader br = null;
		try{
			br = new BufferedReader(new InputStreamReader(new FileInputStream(filename)));
		}catch(Exception e){
			e.printStackTrace();
			return null;
		}
		String line = null;
		Sequence seq1;
		Sequence seq2;
		String seq1tmp = "";
		String seq2tmp = "";
		String ali1tmp;
		String ali2tmp;
		String[] tmp = new String[3];
		try{
			int count = 0;
			while((line = br.readLine()) != null){
				tmp[count%3] = line;
				count++;
				if(count%3 == 0){
					ali1tmp = tmp[1].split(":")[1].trim();
					ali2tmp = tmp[2].split(":")[1].trim();
					for(int i = 0; i != ali1tmp.length(); i++){
						if(ali1tmp.charAt(i) != '-'){
							seq1tmp += ali1tmp.charAt(i);
						}
						if(ali2tmp.charAt(i) != '-'){
							seq2tmp += ali2tmp.charAt(i);
						}
					}
					alignments.add(new SequenceAlignment(new Sequence(tmp[1].split(":")[0].trim(),seq1tmp), new Sequence(tmp[2].split(":")[0].trim(),seq2tmp), tmp[1].split(":")[1].trim(), tmp[2].split(":")[1].trim(), (int)(Gotoh.FACTOR*Double.parseDouble(tmp[0].split("\\s")[2].trim()))));
				}
			}
		} catch(Exception e){
			e.printStackTrace();
		}
		return alignments;
	}
	
	public List<SequenceAlignment> readAlignments(String filename){
		this.filename = filename;
		return readAlignments();
	}
	
	public void initSequentialRead(){
		if(filename != null && !filename.equals("")){
			br = null;
		}
		try{
			br = new BufferedReader(new InputStreamReader(new FileInputStream(filename)));
		}catch(Exception e){
			e.printStackTrace();
		}
	}
	
	public void initSequentialRead(String filename){
		filename = filename;
		try{
			br = new BufferedReader(new InputStreamReader(new FileInputStream(filename)));
		}catch(Exception e){
			e.printStackTrace();
		}
	}
	
	public SequenceAlignment nextAlignment(){
		if(br == null){
			return null;
		}
		String tmp[] = new String[3];
		String line;
		try{
			for(int i = 0; i != 3; i++){
				if((line = br.readLine()) != null){
					tmp[i] = line.trim();
				}else{
					return null;
				}
			}
		}catch(Exception e){
			e.printStackTrace();
			return null;
		}
		
		String ali1tmp = tmp[1].split(":")[1].trim();
		String ali2tmp = tmp[2].split(":")[1].trim();
		String seq1tmp = "";
		String seq2tmp = "";
		for(int i = 0; i != ali1tmp.length(); i++){
			if(ali1tmp.charAt(i) != '-'){
				seq1tmp += ali1tmp.charAt(i);
			}
			if(ali2tmp.charAt(i) != '-'){
				seq2tmp += ali2tmp.charAt(i);
			}
		}
		return new SequenceAlignment(new Sequence(tmp[1].split(":")[0].trim(),seq1tmp), new Sequence(tmp[2].split(":")[0].trim(),seq2tmp), tmp[1].split(":")[1].trim(), tmp[2].split(":")[1].trim(), (int)(Gotoh.FACTOR*Double.parseDouble(tmp[0].split("\\s")[2].trim())));

	}
	
	
}
