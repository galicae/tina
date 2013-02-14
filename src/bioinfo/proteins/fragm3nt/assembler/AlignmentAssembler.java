package bioinfo.proteins.fragm3nt.assembler;

import java.util.LinkedList;

import bioinfo.Sequence;
import bioinfo.alignment.SequenceAlignment;
import bioinfo.alignment.gotoh.FreeshiftSequenceGotoh;
import bioinfo.alignment.matrices.QuasarMatrix;
import bioinfo.proteins.PDBEntry;
import bioinfo.proteins.PDBFileReader;
import bioinfo.proteins.fragm3nt.ProteinFragment;
import bioinfo.superpos.PDBReduce;

public class AlignmentAssembler extends Assembler {

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
		String query = seqs.get(0).getSequenceAsString();
		PDBFileReader reader = new PDBFileReader(pdbDirectory);
		LinkedList<ProteinFragment> result = new LinkedList<ProteinFragment>();

		// read PDBs as ProteinFragments
		LinkedList<ProteinFragment> structures = new LinkedList<ProteinFragment>();
		for (int i = 0; i < seqs.size(); i++) {
			String id = seqs.get(i).getID();
			PDBEntry temp = reader.readFromFolderById(id);
			ProteinFragment tempFrag = new ProteinFragment(temp.getID()
					+ temp.getChainID() + temp.getChainIDNum(),
					PDBReduce.reduceSinglePDB(temp), extent);
			tempFrag.setSequence(query);
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

		ProteinFragment resultFragment = assembleProtein(result, extent, seqs
				.get(0).getSequenceAsString());
		return resultFragment;
	}

	protected LinkedList<ProteinFragment> collectFragmentsFromAl(
			LinkedList<SequenceAlignment> alignments,
			LinkedList<ProteinFragment> structures, int extent) {

		// define query only once - no need to spend three lines every time
		String query = alignments.get(0).getComponent(0).getSequenceAsString();

		// finding curFrag now works with the position and not the sequence
		ProteinFragment curFrag = findFragFromAli(0, alignments, structures,
				null);
		// first find first fragment
		LinkedList<ProteinFragment> result = new LinkedList<ProteinFragment>();
		LinkedList<ProteinFragment> positSolutions = new LinkedList<ProteinFragment>();
		result.add(curFrag);
		int add = fragLength - extent;

		// now loop over bulk of protein
		for (int i = add; i < query.length() - fragLength; i += add) {
			positSolutions.clear();
			curFrag = findFragFromAli(i, alignments, structures, curFrag);

			result.add(curFrag);
		}
		// and now the last
		curFrag = findFragFromAli(query.length() - fragLength, alignments,
				structures, curFrag);
		result.add(curFrag);
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

		double[][] matrix = QuasarMatrix.DAYHOFF_MATRIX;
		String query = "";
		/*
		 * for all alignments, check if the fragment starting (in the query) at
		 * start with length fragLength (no gaps) has a 100% match. If yes, then
		 * score it. If not, save the score as negative infinity.
		 */
		char[][] ali = new char[1][1];
		int[][] map = new int[1][1];
		SequenceAlignment temp = alignments.get(0);
		double[] scores = new double[alignments.size()];
		for (int i = 0; i < alignments.size(); i++) {
			temp = alignments.get(i);
			ali = temp.getRows();
			map = temp.getAlignedResidues();
			int diff1 = map[0][start + fragLength - 1] - map[0][start];
			int diff2 = map[1][start + fragLength - 1] - map[1][start];
			// if one of the two sequences contains a gap continue
			if (diff1 != fragLength || diff2 != fragLength) {
				scores[i] = Double.NEGATIVE_INFINITY;
				continue;
			}
			// else score the aligned part
			for (int j = 0; j < fragLength; j++) {
				char x = ali[0][map[0][start + j]];
				char y = ali[1][map[1][start + j]];
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
		 * if not all scores are negative (a continuous match of length 8 has
		 * been found in the alignments) then return this fragment.
		 */
		if (maxIndex == -1) {
			ProteinFragment result = structures.get(0).getPart(start,
					start + fragLength);
			return result;
		} else {
			map = alignments.get(maxIndex).getAlignedResidues();
			ProteinFragment result = structures.get(maxIndex).getPart(
					map[1][start], map[1][start + fragLength]);
			return result;
		}
	}
}