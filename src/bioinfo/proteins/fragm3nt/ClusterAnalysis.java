package bioinfo.proteins.fragm3nt;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.List;

import bioinfo.proteins.AminoAcid;
import bioinfo.proteins.AminoAcidName;
import bioinfo.proteins.Atom;
import bioinfo.proteins.PDBEntry;

public class ClusterAnalysis {
	private ArrayList<FragmentCluster> clusters;

	public ClusterAnalysis(ArrayList<FragmentCluster> clusters) {
		this.clusters = new ArrayList<FragmentCluster>();
		for (FragmentCluster f : clusters) {
			this.clusters.add(f);
		}
	}

	public ClusterAnalysis(String folder) {
		this.clusters = new ArrayList<FragmentCluster>();
		BufferedReader br = null;
		String id;
		int l = 0;
		try {
			if (folder == null) {
				System.err.println("null folder - no clusters read");
			}
			File[] files = new File(folder).listFiles();
			if (files.length == 0) {
				System.err.println("no files found - no clusters read");
			}
			for (File file : files) {
				br = new BufferedReader(new InputStreamReader(
						new FileInputStream(file)));
				l = file.getName().split("_").length;
				id = file.getName().split("_")[l - 1];
				clusters.add(parseCluster(br, id));
			}
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
	
	private FragmentCluster parseCluster(BufferedReader br, String clusterId) {
		try {
			FragmentCluster result = new FragmentCluster();
			String line = br.readLine();
			String seq = "";
			int fragLength = 0;
			ArrayList<Atom> atoms = new ArrayList<Atom>();
			String name;
			double[] coord = new double[3];
			while(line != null) {
				if(line.startsWith("REMARK")) {
					seq = line.split(" ")[2];
				}
				else if(line.startsWith("MODEL")) {
					clusterId += line.split("\\s+")[1];
				}
				else if(line.startsWith("ATOM")) {
					name = line.substring(12,16).trim();
					coord[0] = Double.parseDouble(line.substring(30,38).trim());
					coord[1] = Double.parseDouble(line.substring(38,46).trim());
					coord[2] = Double.parseDouble(line.substring(46,54).trim());
					atoms.add(new Atom(name,coord));
				}
				else if(line.startsWith("ENDMDL")) {
					Atom[] atomArray = new Atom[atoms.size()];
					for(int i = 0; i < atoms.size(); i++) {
						Atom tempAt = atoms.get(i);
						atomArray[i] = new Atom(tempAt.getType(), tempAt.getPosition());
					}
					ProteinFragment heyIJustMetYou = new ProteinFragment(clusterId, seq, atomArray, fragLength);
					result.add(heyIJustMetYou);
					result.add(heyIJustMetYou);
				}
				
				line = br.readLine();
			}
		}
		catch (Exception e) {
			e.printStackTrace();
		}
		
		return null;
	}

}
