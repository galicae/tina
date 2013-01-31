package test;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.util.LinkedList;
import java.util.List;

import bioinfo.proteins.PDBEntry;
import bioinfo.proteins.PDBFileReader;
import bioinfo.proteins.fragm3nt.Assembler;
import bioinfo.proteins.fragm3nt.FragmentCluster;
import bioinfo.proteins.fragm3nt.Fragmenter;
import bioinfo.proteins.fragm3nt.KMeansAllvsAll;
import bioinfo.proteins.fragm3nt.ProteinFragment;
import bioinfo.superpos.Kabsch;
import bioinfo.superpos.PDBReduce;
import bioinfo.superpos.Transformation;

public class Fr4gmentTest {
	static int fragLength = 6;

	public static void main(String[] args) {
		PDBFileReader reader = new PDBFileReader("./proteins/");
		List<PDBEntry> list = new LinkedList<PDBEntry>();
		list = reader.readPdbFolder();
		double cutoff = 0.01;
		PDBEntry pdb1 = list.get(0);
		
		for(int i = 0; i < list.size(); i++) {
			pdb1 = list.get(i);
			if(hasNoRepeats(pdb1))
				magic(pdb1, cutoff);
		}
	}

	
	public static boolean hasNoRepeats(PDBEntry pdb) {
		LinkedList<String> frags = new LinkedList<String>();
		String curFrag = "";
		for(int i = 0; i < pdb.length() - fragLength; i++) {
			frags.add("");
			for (int j = i; j < i + fragLength; j++) {
				curFrag += pdb.getAminoAcid(j).getName().getOneLetterCode();
			}
			// check if frag is already there
			for(int k = 0; k < frags.size(); k++) {
				if(curFrag.equals(frags.get(k)))
					return false;
			}
			frags.add(curFrag);
			curFrag = "";
		}
		return true;
	}
	
	
	public static void magic(PDBEntry pdb1, double cutoff) {
		System.out.println("Magicking protein " + pdb1.getID());
		LinkedList<ProteinFragment> pList = new LinkedList<ProteinFragment>();
		Fragmenter.crunchBackboneSeq(pdb1, pList, fragLength);
		KMeansAllvsAll clustah = new KMeansAllvsAll(pList, cutoff);
		LinkedList<FragmentCluster> clusters = new LinkedList<FragmentCluster>();
		clustah.initializeClusters(0);
		clustah.update(0);

		clusters = clustah.getClusters();

		// assembly?
		Assembler ass = new Assembler(fragLength);
		String query = ass.readSequence(pdb1);
		ProteinFragment resultFragment = ass.predictStructure(query, clusters, 3);

		double[][][] kabschFood = new double[2][pdb1.length()][3];
		kabschFood[0] = PDBReduce.reduceSinglePDB(pdb1);

		ProteinFragment prot = new ProteinFragment("real", query,
				new double[0][0], fragLength);
		prot.append(kabschFood[0], "");

		kabschFood[1] = resultFragment.getAllResidues();
		Transformation t = Kabsch.calculateTransformation(kabschFood);

		resultFragment.setCoordinates(t.transform(resultFragment
				.getAllResidues()));
		System.out.println(t.getRmsd());

		try {
			BufferedWriter w = new BufferedWriter(new FileWriter("./results/sh"
					+ pdb1.getID() + ".pdb"));
			w.write("MODEL        1\n");
			w.write(resultFragment.toString());
			w.write("ENDMDL\n");
			w.write("MODEL        2\n");
			w.write(prot.toString());
			w.write("ENDMDL\n");
			w.close();
		} catch (Exception e) {
			e.printStackTrace();
		}
	}

}
