package bioinfo.proteins.fragm3nt.assembler;

import java.util.LinkedList;

import bioinfo.Sequence;
import bioinfo.alignment.SequenceAlignment;
import bioinfo.alignment.gotoh.FreeshiftSequenceGotoh;
import bioinfo.alignment.matrices.QuasarMatrix;
import bioinfo.proteins.PDBEntry;
import bioinfo.proteins.PDBFileReader;
import bioinfo.proteins.fragm3nt.ProteinFragment;
import bioinfo.superpos.Kabsch;
import bioinfo.superpos.PDBReduce;
import bioinfo.superpos.Transformation;

public class AlignmentAssembler extends Assembler {
	private LinkedList<ProteinFragment> frags = new LinkedList<ProteinFragment>();

	public AlignmentAssembler(int fragLength) {
		super(fragLength);
	}

	/**
	 * 
	 * @param query
	 * @param seqs
	 *            convention: first sequence is the query
	 * @param extent
	 * @return
	 */
	public ProteinFragment predStrucFromAl(LinkedList<Sequence> seqs,
			int extent, String pdbDirectory) {

		PDBFileReader reader = new PDBFileReader(pdbDirectory);
		LinkedList<ProteinFragment> result = new LinkedList<ProteinFragment>();

		// read PDBs as ProteinFragments
		LinkedList<ProteinFragment> structures = new LinkedList<ProteinFragment>();
		for (int i = 0; i < seqs.size(); i++) {
			String id = seqs.get(i).getID();
			// query = seqs.get(i).getSequenceAsString();
			PDBEntry temp = reader.readFromFolderById(id);
			ProteinFragment tempFrag = new ProteinFragment(temp.getID(),
					PDBReduce.reduceSinglePDB(temp), extent);
			tempFrag.setSequence(temp.getSequence());
			structures.add(tempFrag);
		}

		// create alignments
		LinkedList<SequenceAlignment> alignments = new LinkedList<SequenceAlignment>();
		FreeshiftSequenceGotoh fg = new FreeshiftSequenceGotoh(-15, -3,
				QuasarMatrix.DAYHOFF_MATRIX);
		for (int i = 1; i < seqs.size(); i++) {
			alignments.add(fg.align(seqs.get(0), seqs.get(i)));
		}

		// collect fragments
		result = collectFragmentsFromAl(alignments, structures, extent);
		for (int i = 0; i < result.size(); i++) {
			frags.add(result.get(i).clone());
		}
		ProteinFragment tempResFragment = assembleProtein(result, extent, seqs
				.get(0).getSequenceAsString());
		ProteinFragment resultFragment = new ProteinFragment("RE",
				new double[1][1], extent);
		resultFragment = tempResFragment.clone();
		resultFragment.setID(seqs.get(0).getID());
		resultFragment.setClusterIndex(result.get(0).getClusterIndex());
		return resultFragment;
	}

	/**
	 * 
	 * @param alignments
	 * @param structures
	 * @param extent
	 * @return
	 */
	protected LinkedList<ProteinFragment> collectFragmentsFromAl(
			LinkedList<SequenceAlignment> alignments,
			LinkedList<ProteinFragment> structures, int extent) {
		int count = 0;
//		System.out.println(alignments.get(0).toStringVerbose());
		// define query only once - no need to spend three lines every time
		String query = alignments.get(0).getComponent(0).getSequenceAsString();

		// finding curFrag now works with the position and not the sequence
		ProteinFragment curFrag = findFragFromAli(0, alignments, structures,
				null);
		// first find first fragment
		LinkedList<ProteinFragment> result = new LinkedList<ProteinFragment>();
		LinkedList<ProteinFragment> positSolutions = new LinkedList<ProteinFragment>();
		result.add(curFrag);
		if (!curFrag.getID().startsWith(structures.get(0).getID())) {
			count++;
		}
		int add = fragLength - extent;

		// now loop over bulk of protein
		for (int i = add; i < query.length() - fragLength;) {
			positSolutions.clear();
			curFrag = findFragFromAli(i, alignments, structures, curFrag);
			if (!curFrag.getID().startsWith(structures.get(0).getID())) {
				count++;
			}
			result.add(curFrag);
			i += add;
		}
		// and now the last
		curFrag = findFragFromAli(query.length() - fragLength - 1, alignments,
				structures, curFrag);
		result.add(curFrag);
		if (!curFrag.getID().startsWith(structures.get(0).getID())) {
			count++;
		}
		result.get(0).setClusterIndex((int)(count / (1.0 * result.size()) * 1000));
		return result;
	}

	/**
	 * 
	 * @param start
	 * @param alignments
	 * @param structures
	 * @param lastFrag
	 * @return
	 */
	protected ProteinFragment findFragFromAli(int start,
			LinkedList<SequenceAlignment> alignments,
			LinkedList<ProteinFragment> structures, ProteinFragment lastFrag) {

		int saveStart = start;
		double[][] matrix = QuasarMatrix.DAYHOFF_MATRIX;
		/*
		 * for all alignments, check if the fragment starting (in the query) at
		 * start with length fragLength (no gaps) has a 100% match. If yes, then
		 * score it. If not, save the score as negative infinity.
		 */
		int[][] map = new int[1][1];
		SequenceAlignment temp = alignments.get(0);
		double[] scores = new double[alignments.size()];
		for (int i = 0; i < alignments.size(); i++) {

			temp = alignments.get(i);
			map = temp.getAlignedResidues();
			// find start:
			boolean found = false;
			for (int j = 0; j < map[0].length; j++) {
				if (map[0][j] == start) {
					start = j;
					found = true;
					break;
				}
			}
			if (start < map[0].length - fragLength - 1 && found) {
				int diff1 = map[0][start + fragLength - 1] - map[0][start] + 1;
				int diff2 = map[1][start + fragLength - 1] - map[1][start] + 1;

				// if one of the two sequences contains a gap continue
				if (diff1 != fragLength || diff2 != fragLength) {
					scores[i] = Double.NEGATIVE_INFINITY;
					continue;
				}
			} else {
				scores[i] = Double.NEGATIVE_INFINITY;
				continue;
			}
			// else score the aligned part
			for (int j = 0; j < fragLength; j++) {
				char x = temp.getComponent(0).getSequence()[map[0][start + j]];
				char y = temp.getComponent(1).getSequence()[map[1][start + j]];
				scores[i] += matrix[x - 65][y - 65];
			}
		}

		/*
		 * if all partial scores are all negative infinity, just copy the
		 * fragment from the native.
		 */
		double maxScore = Double.NEGATIVE_INFINITY;
		int maxIndex = -1;
		for (int i = 0; i < scores.length; i++) {
			if (scores[i] > maxScore) {
				maxIndex = i;
				maxScore = scores[i];
			}
		}

		/*
		 * if not all scores are negative (a continuous match of length fragLength has
		 * been found in the alignments) then return this fragment.
		 */
		if (maxIndex == -1) {
			ProteinFragment result = structures.get(0).getPart(saveStart,
					saveStart + fragLength);
//			System.out.println(result.getID() + " " + result.getSequence());
			return result;
		} else {
			map = alignments.get(maxIndex).getAlignedResidues();
			ProteinFragment result = structures.get(maxIndex + 1).getPart(
					map[1][start], map[1][start] + fragLength);
			if(checkFragment(result)) {
//			if(true) {
//				System.out.println(result.getID() + " " + result.getSequence());
				return result;
			}
			else
				result = structures.get(0).getPart(saveStart,
						saveStart + fragLength);
//			System.out.println(result.getID() + " " + result.getSequence());
			return result;
		}
	}

	protected boolean checkFragment(ProteinFragment f) {
		double dist = 0;
		for (int i = 1; i < f.fragLength; i++) {
			dist = distance(f.getResidue(i - 1), f.getResidue(i));
			if (dist < 3.5 || dist > 4.0) {
				return false;
			}
		}
		return true;
	}

	protected double distance(double[] a, double[] b) {
		double sum = 0;
		for (int i = 0; i < a.length; i++) {
			sum += (a[i] - b[i]) * (a[i] - b[i]);
		}
		return Math.sqrt(sum);
	}

	protected void alignFragments(ProteinFragment stable, ProteinFragment move,
			int extent) {
		double[][][] kabschFood = new double[2][extent][3];
		int shove = stable.getAllResidues().length - extent;
		for (int i = 0; i < extent; i++) {
			kabschFood[0][i] = stable.getResidue(shove + i);
			kabschFood[1][i] = move.getResidue(i);
		}
		Transformation t = Kabsch.calculateTransformation(kabschFood);
		move.setCoordinates(t.transform(move.getAllResidues()));
		matchCoordinates(stable, move, extent);

		double[][] moveCoordinates = new double[(fragLength - extent)][3];
		double[][] allResidues = move.getAllResidues();
		for (int i = 0; i < moveCoordinates.length; i++) {
			moveCoordinates[i][0] = allResidues[extent + i][0];
			moveCoordinates[i][1] = allResidues[extent + i][1];
			moveCoordinates[i][2] = allResidues[extent + i][2];
		}
		String moveSeq = move.getSequence().substring(extent);
		stable.append(moveCoordinates, moveSeq);
		// System.out.println();
	}

	public LinkedList<ProteinFragment> getFragments() {
		return this.frags;
	}
}
