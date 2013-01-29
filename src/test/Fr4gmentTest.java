package test;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.util.LinkedList;
import java.util.List;

import bioinfo.proteins.Atom;
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
	static int fragLength = 5;
	public static void main(String[] args) {
		PDBFileReader reader = new PDBFileReader("./proteins2/");
		List<PDBEntry> files = new LinkedList<PDBEntry>();
		LinkedList<ProteinFragment> pList = new LinkedList<ProteinFragment>();
		files = reader.readPdbFolder();
		// files.add(pdb1);
		PDBEntry pdb1 = files.get(0);
		for (PDBEntry e : files) {
			Fragmenter.crunchBackboneSeq(e, pList, fragLength);
		}
		System.out.println("crunched");
		KMeansAllvsAll clustah = new KMeansAllvsAll(pList, 0.01);
		LinkedList<FragmentCluster> clusters = new LinkedList<FragmentCluster>();
		clustah.initializeClusters();
		clustah.update(0);

		clusters = clustah.getClusters();
		
		// assembly?
		Assembler ass = new Assembler(fragLength);
		String query = ass.readSequence(pdb1);
		ProteinFragment resultFragment = ass.predictStructure(query, clusters, 3);
		
		
		double[][][] kabschFood = new double[2][pdb1.length()][3];
		kabschFood[0] = PDBReduce.reduceSinglePDB(pdb1);
		
		ProteinFragment prot = new ProteinFragment("real", query, new double[0][0], fragLength);
		prot.append(kabschFood[0], "");
		
		kabschFood[1] = resultFragment.getAllResidues();
		Transformation t = Kabsch.calculateTransformation(kabschFood);
		
		resultFragment.setCoordinates(t.transform(resultFragment.getAllResidues()));
		System.out.println(t.getRmsd());
		
		try {
			BufferedWriter w = new BufferedWriter(new FileWriter("C:/users/nikos/Desktop/test.pdb"));
			w.write("MODEL        1\n");
			w.write(resultFragment.toString());
			w.write("ENDMDL\n");
			w.write("MODEL        2\n");
			w.write(prot.toString());
			w.write("ENDMDL\n");
			w.close();
		}
		catch(Exception e) {
			e.printStackTrace();
		}
	}

}
