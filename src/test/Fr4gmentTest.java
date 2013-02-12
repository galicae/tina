package test;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.util.LinkedList;
import java.util.List;

import bioinfo.proteins.PDBEntry;
import bioinfo.proteins.PDBFileReader;
import bioinfo.proteins.fragm3nt.ClusterAnalysis;
import bioinfo.proteins.fragm3nt.FragmentCluster;
import bioinfo.proteins.fragm3nt.*;
import bioinfo.superpos.Kabsch;
import bioinfo.superpos.PDBReduce;
import bioinfo.superpos.Transformation;

public class Fr4gmentTest {

	public static void main(String[] args) {
		PDBFileReader reader = new PDBFileReader("./proteins2/");
		List<PDBEntry> fileList = new LinkedList<PDBEntry>();
		fileList = reader.readPdbFolder();
		double sum = 0;
		int fragLength = 8;
		int extent = 5;
		String query = "";

		ClusterAnalysis c = new ClusterAnalysis("clusters");
		LinkedList<FragmentCluster> clusters = c.getClusters();

		double[][][] kabschFood;

		try {
			BufferedWriter wr = new BufferedWriter(new FileWriter(
					"famTestResults"));
			for (PDBEntry pdb : fileList) {
				CheatAssembler ass = new CheatAssembler(fragLength, pdb);
				BufferedWriter wr2 = new BufferedWriter(new FileWriter(
						"famTest/" + pdb.getID() + ".pdb"));
				query = "";
				for (int i = 0; i < pdb.length(); i++) {
					query += pdb.getAminoAcid(i).getName().getOneLetterCode();
				}
				try {
					ProteinFragment result = ass.proveConcept(query, clusters,
							extent);
					kabschFood = new double[2][result.getAllResidues().length][3];
					kabschFood[0] = result.getAllResidues();
					kabschFood[1] = PDBReduce.reduceSinglePDB(pdb);
					if (kabschFood[1] == null)
						continue;
					Transformation t = Kabsch
							.calculateTransformation(kabschFood);
					result.setCoordinates(t.transform(kabschFood[0]));
					wr.write(fileList.indexOf(pdb) + "\t" + pdb.getID() + "\t" + Double.toString(t.getRmsd()) + "\n");
					System.out.println(t.getRmsd());
					sum += t.getRmsd();

					ProteinFragment prot = new ProteinFragment("real", query, new double[0][0], fragLength);
					prot.append(PDBReduce.reduceSinglePDB(pdb), "");
					prot.setCoordinates(kabschFood[1]);

					wr2.write("MODEL        1\n");
					wr2.write(result.toString());
					wr2.write("ENDMDL\n");
					wr2.write("MODEL        2\n");
					wr2.write(prot.toString());
					wr2.write("ENDMDL\n");
					wr2.close();
				} catch (Exception e) {
					continue;
				}
			}
			System.out.println(sum / fileList.size());
			wr.close();
		} catch (Exception e) {
			e.printStackTrace();
		}

	}
	// public static boolean hasNoRepeats(PDBEntry pdb) {
	// LinkedList<String> frags = new LinkedList<String>();
	// String curFrag = "";
	// for (int i = 0; i < pdb.length() - fragLength; i++) {
	// frags.add("");
	// for (int j = i; j < i + fragLength; j++) {
	// curFrag += pdb.getAminoAcid(j).getName().getOneLetterCode();
	// }
	// // check if frag is already there
	// for (int k = 0; k < frags.size(); k++) {
	// if (curFrag.equals(frags.get(k)))
	// return false;
	// }
	// frags.add(curFrag);
	// curFrag = "";
	// }
	// return true;
	// }
	//
	// public static void magic(PDBEntry pdb1, double cutoff) {
	// System.out.println("Magicking protein " + pdb1.getID());
	// LinkedList<ProteinFragment> pList = new LinkedList<ProteinFragment>();
	// Fragmenter.crunchBackboneSeq(pdb1, pList, fragLength);
	// KMeansAllvsAll clustah = new KMeansAllvsAll(pList, cutoff);
	// LinkedList<FragmentCluster> clusters = new LinkedList<FragmentCluster>();
	// clustah.initializeClusters(1);
	// clustah.update(10);
	//
	// clusters = clustah.getClusters();
	//
	// // assembly?
	// Assembler ass = new Assembler(fragLength);
	// String query = ass.readSequence(pdb1);
	// ProteinFragment resultFragment = ass.predictStructure(query, clusters,
	// 4);
	//
	//
	// double[][][] kabschFood = new double[2][pdb1.length()][3];
	// kabschFood[0] = PDBReduce.reduceSinglePDB(pdb1);
	// ProteinFragment prot = new ProteinFragment("real", query,
	// new double[0][0], fragLength);
	//
	// if(kabschFood[0] == null) {
	// System.out.println("couldn't magic protein " + pdb1.getID());
	// return;
	// }
	//
	// prot.append(PDBReduce.reduceSinglePDB(pdb1), "");
	//
	// kabschFood[1] = resultFragment.getAllResidues();
	// Transformation t = Kabsch.calculateTransformation(kabschFood);
	//
	// resultFragment.setCoordinates(t.transform(resultFragment
	// .getAllResidues()));
	// rmsdResults.add(t.getRmsd());
	//
	// try {
	// BufferedWriter w = new BufferedWriter(new FileWriter("./results/TEST"
	// + pdb1.getID() + ".pdb"));
	// w.write("MODEL        1\n");
	// w.write(resultFragment.toString());
	// w.write("ENDMDL\n");
	// w.write("MODEL        2\n");
	// w.write(prot.toString());
	// w.write("ENDMDL\n");
	// w.close();
	// } catch (Exception e) {
	// e.printStackTrace();
	// }
	// }

}
