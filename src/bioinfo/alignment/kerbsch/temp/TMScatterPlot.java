package bioinfo.alignment.kerbsch.temp;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;

import bioinfo.Sequence;
import bioinfo.alignment.Aligner;
import bioinfo.alignment.SequenceAlignment;
import bioinfo.alignment.gotoh.FreeshiftSequenceGotoh;
import bioinfo.alignment.gotoh.GlobalSequenceGotoh;
import bioinfo.alignment.gotoh.LocalSequenceGotoh;
import bioinfo.alignment.matrices.QuasarMatrix;
import bioinfo.proteins.PDBEntry;
import bioinfo.proteins.PDBFileReader;
import bioinfo.superpos.TMMain;
import bioinfo.superpos.Transformation;

/*
 * This class creates a scatter-plot for 2 alignment methods ranked by their TM-score
 */

public class TMScatterPlot {
	private TMMain superpos = new TMMain();
	private HashMap<String, char[]> seqlib1;
	private HashMap<String, char[]> seqlib2;
	private ArrayList<String[]> pairs;
	private Aligner method1;
	private Aligner method2;
	private BufferedWriter out;
	private PDBFileReader pdbreader;

	public TMScatterPlot(Aligner method1, Aligner method2,
			HashMap<String, char[]> seqlib1, HashMap<String, char[]> seqlib2, ArrayList<String[]> pairs,
			String pdbpath, String outpath) {
		this.method1 = method1;
		this.method2 = method2;
		this.seqlib1 = seqlib1;
		this.seqlib2 = seqlib2;
		this.pairs = pairs;
		this.pdbreader = new PDBFileReader(pdbpath);
		try {
			this.out = new BufferedWriter(new FileWriter(outpath));
		} catch (IOException e) {
			System.out.println("cannot initialize writer(scatterplot)");
		}
	}

	private void writeTMScores(String query, String template,
			SequenceAlignment align1, SequenceAlignment align2)
			throws Exception {
		PDBEntry p = pdbreader.readFromFolderById(query);
		PDBEntry q = pdbreader.readFromFolderById(template);
		Transformation tr1 = superpos.calculateTransformation(align1, p, q);
		Transformation tr2 = superpos.calculateTransformation(align2, p, q);
		out.append(tr1.getTmscore() + "\t" + tr2.getTmscore() + "\n");
	}

	public void makePlot() {
		Sequence query1;
		Sequence template1;
		Sequence query2;
		Sequence template2;
		SequenceAlignment align1;
		SequenceAlignment align2;
		
		for (String[] pair : this.pairs) {
			query1 = new Sequence(pair[0],seqlib1.get(pair[0]));
			template1 = new Sequence(pair[1],seqlib1.get(pair[1]));
			query2 = new Sequence(pair[0],seqlib2.get(pair[0]));
			template2 = new Sequence(pair[1],seqlib2.get(pair[1]));
			align1 = (SequenceAlignment) method1.align(query1,template1);
			align2 = (SequenceAlignment) method2.align(query2,template2);

			try {
				writeTMScores(pair[0], pair[1], align1, align2);
			} catch (Exception e) {
				System.out
						.println("something is wrong with TM-score calculation in scatterplot");
				e.printStackTrace();
			}
			System.out.println("created point in plot");
		}
		try {
			out.close();
		} catch (IOException e) {
			System.out.println("cannot close writer(scatterplot)");
		}
	}
	
	public static void main(String[] args){
		Aligner method1 = null;
		Aligner method2 = null;
		
		InitClass init = new InitClass();
		//method1
		String feature1 = args[0];
		double go1 = Double.parseDouble(args[1]);
		double ge1 = Double.parseDouble(args[2]);
		String mode1 = args[3];
		double[][] matrix1 = null;
		HashMap<String,char[]> seqlib1 = SeqLibrary.read(args[4]);
		
		//method2
		String feature2 = args[5];
		double go2 = Double.parseDouble(args[6]);
		double ge2 = Double.parseDouble(args[7]);
		String mode2 = args[8];
		double[][] matrix2 = null;
		HashMap<String,char[]> seqlib2 = SeqLibrary.read(args[9]);
		
		//other stuff
		ArrayList<String[]> pairs = PairReader.parse(args[10]);
		String pdbpath = args[11];
		String outpath = args[12];
		
		//create matrices
		if(feature1.equals("polarity")){
			matrix1 = init.calcGotohInputMatrix(init.calcPolarityScores());
		} else if(feature1.equals("hydrophob")){
			matrix1 = init.calcGotohInputMatrix(init.calcHydropathyScores());
		} else if(feature1.equals("secstruct")){
			matrix1 = SecStructScores.matrix;
		} else if(feature1.equals("sequence")){
			matrix1 = QuasarMatrix.DAYHOFF_MATRIX;
		}
		
		if(feature2.equals("polarity")){
			matrix2 = init.calcGotohInputMatrix(init.calcPolarityScores());
		} else if(feature2.equals("hydrophob")){
			matrix2 = init.calcGotohInputMatrix(init.calcHydropathyScores());
		} else if(feature2.equals("secstruct")){
			matrix2 = SecStructScores.matrix;
		} else if(feature2.equals("sequence")){
			matrix2 = QuasarMatrix.DAYHOFF_MATRIX;
		}
		
		//create Aligner
		if(mode1.equals("local") && !feature1.equals("angle")){
			method1 = new LocalSequenceGotoh(go1,ge1,matrix1);
		} else if(mode1.equals("global") && !feature1.equals("angle")){
			method1 = new GlobalSequenceGotoh(go1,ge1,matrix1);
		} else if(mode1.equals("freeshift") && !feature1.equals("angle")){
			method1 = new FreeshiftSequenceGotoh(go1,ge1,matrix1);
		}
		
		if(mode2.equals("local") && !feature2.equals("angle")){
			method2 = new LocalSequenceGotoh(go2,ge2,matrix2);
		} else if(mode2.equals("global") && !feature2.equals("angle")){
			method2 = new GlobalSequenceGotoh(go2,ge2,matrix2);
		} else if(mode2.equals("freeshift") && !feature2.equals("angle")){
			method2 = new FreeshiftSequenceGotoh(go2,ge2,matrix2);
		}
		
		//calculate the scatterplot
		TMScatterPlot plot = new TMScatterPlot(method1,method2,seqlib1,seqlib2,pairs,pdbpath,outpath);
		plot.makePlot();
	}
}
