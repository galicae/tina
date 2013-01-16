package bioinfo.alignment;

import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.List;

import bioinfo.Sequence;
import bioinfo.alignment.gotoh.Gotoh;

/**
 * reads PDBFile and returns it as internal Alignment
 * 
 * @author andreseitz
 * 
 */
public class SequenceAlignmentFileReader {

	private String filename = null;
	private BufferedReader br;

	public SequenceAlignmentFileReader() {

	}

	public SequenceAlignmentFileReader(String filename) {
		this.filename = filename;
	}

	/**
	 * reads alignments from file @filename and returns them in a list
	 * 
	 * @return the list of all read alignments
	 */
	public List<SequenceAlignment> readAlignments() {
		List<int[]> map = new ArrayList<int[]>();
		List<SequenceAlignment> alignments = new ArrayList<SequenceAlignment>();
		BufferedReader br = null;
		try {
			br = new BufferedReader(new InputStreamReader(new FileInputStream(
					filename)));
		} catch (Exception e) {
			e.printStackTrace();
			return null;
		}
		String line = null;
		String seq1tmp = "";
		String seq2tmp = "";
		String ali1tmp;
		String ali2tmp;
		String[] tmp = new String[3];
		try {
			int count = 0;
			while ((line = br.readLine()) != null) {
				tmp[count % 3] = line;
				count++;
				if (count % 3 == 0) {
					ali1tmp = tmp[1].split(":")[1].trim();
					ali2tmp = tmp[2].split(":")[1].trim();
					for (int i = 0; i != ali1tmp.length(); i++) {
						if (ali1tmp.charAt(i) != '-') {
							seq1tmp += ali1tmp.charAt(i);
						}
						if (ali2tmp.charAt(i) != '-') {
							seq2tmp += ali2tmp.charAt(i);
						}
					}
					
					//calculate map
					int x = 0, y = 0;
					for (int i = 0; i < ali1tmp.length(); i++) {
						if(ali1tmp.charAt(i) != '-'){
							x++;
							if(ali2tmp.charAt(i) != '-'){
								y++;
								map.add(new int[]{x,y}); //store aligned indices of the two sequences
								continue;
							}
						}
						if(ali2tmp.charAt(i) != '-'){
							y++;
							if(ali1tmp.charAt(i) != '-'){
								x++;
								map.add(new int[]{x,y}); //store aligned indices of the two sequences
								continue;
							}
						}
					}
					alignments.add(new SequenceAlignment(new Sequence(tmp[1]
							.split(":")[0].trim(), seq1tmp), new Sequence(
							tmp[2].split(":")[0].trim(), seq2tmp), tmp[1]
							.split(":")[1].trim(), tmp[2].split(":")[1].trim(),
							(int) (Gotoh.FACTOR * Double.parseDouble(tmp[0]
									.split("\\s")[2].trim())), map.toArray(new int[map.size()][])));
				}
			}
			br.close();
		} catch (Exception e) {
			e.printStackTrace();
		}
		return alignments;
	}

	/**
	 * reads a file containing many alignments and returns them as an
	 * {@link SequenceAlignment} List. Written to read files in the Gotoh output
	 * format
	 * 
	 * @param filename
	 *            the name of the file to read. Should be a path if the file is
	 *            not on the same level as Java
	 * @return a list of {@link SequenceAlignment} objects
	 */
	public List<SequenceAlignment> readAlignments(String filename) {
		this.filename = filename;
		return readAlignments();
	}

	public void initSequentialRead() {
		if (filename != null && !filename.equals("")) {
			br = null;
		}
		try {
			br = new BufferedReader(new InputStreamReader(new FileInputStream(
					filename)));
		} catch (Exception e) {
			e.printStackTrace();
		}
	}

	/**
	 * initialization of the parameters of the class. Needed when an alignment
	 * file is being read sequentially and not in one go. nextAlignment() should
	 * be called after this.
	 * 
	 * @param filename
	 *            the name of the file to read. Should be a path if the file is
	 *            not on the same level as Java
	 */
	public void initSequentialRead(String filename) {
		this.filename = filename;
		try {
			br = new BufferedReader(new InputStreamReader(new FileInputStream(
					filename)));
		} catch (Exception e) {
			e.printStackTrace();
		}
	}

	/**
	 * retrieves the next alignment out of a sequence alignment file; called if
	 * the alignment file would create a list liable to fill the memory and
	 * considerably slow processes. Makes lots of sense when the alignment is
	 * somehow processed (intensively) and memory should be free for these
	 * processes.
	 * 
	 * @return the next alignment; this object essentially contains the
	 *         sequences without gaps, since they are meant to be used in a
	 *         realignment procedure
	 */
	public SequenceAlignment nextAlignment() {
		List<int[]> map = new ArrayList<int[]>();
		
		if (br == null) {
			return null;
		}
		String tmp[] = new String[3];
		String line;
		try {
			for (int i = 0; i != 3; i++) {
				if ((line = br.readLine()) != null) {
					tmp[i] = line.trim();
				} else {
					return null;
				}
			}
		} catch (Exception e) {
			e.printStackTrace();
			return null;
		}

		String ali1tmp = tmp[1].split(":")[1].trim();
		String ali2tmp = tmp[2].split(":")[1].trim();
		String seq1tmp = "";
		String seq2tmp = "";
		for (int i = 0; i != ali1tmp.length(); i++) {
			if (ali1tmp.charAt(i) != '-') {
				seq1tmp += ali1tmp.charAt(i);
			}
			if (ali2tmp.charAt(i) != '-') {
				seq2tmp += ali2tmp.charAt(i);
			}
		}
		
		//calculate map
		int x = 0, y = 0;
		for (int i = 0; i < ali1tmp.length(); i++) {
			if(ali1tmp.charAt(i) != '-'){
				x++;
				if(ali2tmp.charAt(i) != '-'){
					y++;
					map.add(new int[]{x,y}); //store aligned indices of the two sequences
					continue;
				}
			}
			if(ali2tmp.charAt(i) != '-'){
				y++;
				if(ali1tmp.charAt(i) != '-'){
					x++;
					map.add(new int[]{x,y}); //store aligned indices of the two sequences
					continue;
				}
			}
		}
		return new SequenceAlignment(new Sequence(tmp[1].split(":")[0].trim(),
				seq1tmp), new Sequence(tmp[2].split(":")[0].trim(), seq2tmp),
				tmp[1].split(":")[1].trim(), tmp[2].split(":")[1].trim(),
				(int) (Gotoh.FACTOR * Double.parseDouble(tmp[0].split("\\s")[2]
						.trim())), map.toArray(new int[map.size()][]));

	}

}
