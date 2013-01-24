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

public class Fr4gmentTest {
	static int fragLength = 5;
	public static void main(String[] args) {
		PDBFileReader reader = new PDBFileReader("./proteins2/");
		List<PDBEntry> files = new LinkedList<PDBEntry>();
		LinkedList<ProteinFragment> pList = new LinkedList<ProteinFragment>();
		files = reader.readPdbFolder();
		// files.add(pdb1);
		PDBEntry pdb1 = files.get(1);
		for (PDBEntry e : files) {
			Fragmenter.crunchBackboneSeq(e, pList, fragLength);
		}
		System.out.println("crunched");
		KMeansAllvsAll clustah = new KMeansAllvsAll(pList, 1);
		LinkedList<FragmentCluster> clusters = new LinkedList<FragmentCluster>();
		clustah.initializeClusters();
		clustah.update(20);

		clusters = clustah.getClusters();
		
		// assembly?
		Assembler ass = new Assembler(fragLength);
		String query = ass.readSequence(pdb1);
		ProteinFragment resultFragment = ass.predictStructure(query, clusters, 2);
		
		try {
			BufferedWriter w = new BufferedWriter(new FileWriter("/home/p/papadopoulos/Desktop/test.pdb"));
			w.write(resultFragment.toString());
			w.close();
		}
		catch(Exception e) {
			e.printStackTrace();
		}
	}

}
