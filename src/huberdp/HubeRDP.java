/******************************************************************************
 * huberdp.HubeRDP.java                                                       *
 * This file contains the class HubeRDP which is an RDP implementation for    *
 * the GoBi 2012/13.                                                          *
 *                                                                            *
 * This file is best read at line width 80 and tab width 4.                   *
 *                                                                   huberste *
 ******************************************************************************/
package huberdp;

import bioinfo.Sequence;
import bioinfo.alignment.SequenceAlignment;
import huberdp.oracles.TinyOracle;
import huberdp.scoring.SimpleScoring;

import java.util.LinkedList;

/**
 * @author huberste
 * @lastchange 2013-02-18
 */
public class HubeRDP {

	private final static int M = 3; // number of different solutions an oracle
									// shall give
									// private final static int N = 2; // number
									// of subproblems after using oracle

	private LinkedList<Oracle> oracles;

	private Scoring scoring;

	/**
	 * constructs an HubeRDP
	 */
	public HubeRDP() {
		oracles = new LinkedList<Oracle>();
		scoring = new SimpleScoring();
	}

	/**
	 * adds an oracle to the List of oracles
	 * 
	 * @param oracle
	 */
	public void addOracle(Oracle oracle) {
		oracles.add(oracle);
	}

	/**
	 * first rdp must be called with t = new RDPSolutionTree(); pq = new
	 * RDPPriorityQueue(t.getRoot()); rdp (t, pq); Optimal solution is now in
	 * t.getRoot();
	 */
	// RDP (T, pq):=
	public void rdp(RDPSolutionTree t, RDPPriorityQueue pq) {
		// if (pq = {} ) do return root
		if (pq.isEmpty()) {
			// nothing to do here: everything is already calculated
			// DO NOTHING
		} else {
			// v := <SP, \empty >^{\vee} \leftarrow first(pq)
			RDPSolutionTreeOrNode v = pq.getFirst();
			// U := {<PA, \empty >^{\wedge}}
			RDPSolutionTreeAndNode[] uSet = gAND(v, M, t);
			// U <-- sf_{\wedge}(U)
			sf(uSet);
			// if (U = {}) do <SP, TA>^{\wedge} <-- Finish(v,T)
			if (uSet.length == 0) {
				finish(v, t);
			} else {
				// T <-- insert (T,U)
				v.addChildren(uSet);
				// for each u := <PA, {}>^{\wedge} \in U do
				for (RDPSolutionTreeAndNode u : uSet) {
					// if (Leaf(u)) do <PA,TA>^{\wedge} <-- Finish (u, T)
					// TODO this would always be correct, something's wrong here
					// if (u.isLeaf()) {
					// finish(u,t);
					// } else {
					// V:= {<SP', {}>^{\vee}} <-- g_{\vee}(u, T)
					RDPSolutionTreeOrNode[] vSet = gOR(u, t);

					if (vSet.length == 0) {
						finish(u, t);
					} else {
						// T <-- insert (T, V)
						u.addChildren(vSet);
						// pq <-- insert(pq, V)
						pq.add(vSet);
					}

				}

			}
			rdp(t, pq);
		}
	}

	/**
	 * Function that calls the oracles:<br />
	 * "the generating function g_{\wedge}: V_{\vee} \rightarrow V^m_{\wedge} is
	 * applied to a Or-node v \in V_{\vee} and produces a set of m alternate
	 * nodes u_i \in v_{\wedge}, which represent possible alternative
	 * extensions/completions of the current (partial) alignment coded by v."<br />
	 * (From: Protein Threading by Recursive Dynamic Programming. JMB 290,
	 * 757-779)
	 * 
	 * @param v
	 *            the OrNode of the SolutionTree to be worked on
	 * @param m
	 *            maximum number of alternative extensions this function shall
	 *            create
	 * @param t
	 *            the complete SolutionTree
	 * @return (partial) solutions for this subproblem
	 */
	private RDPSolutionTreeAndNode[] gAND(RDPSolutionTreeOrNode v, int m,
			RDPSolutionTree t) {

		LinkedList<RDPSolutionTreeAndNode> results = new LinkedList<RDPSolutionTreeAndNode>();

		// BEGIN DEBUGGING
		// System.out.println("sending problem to oracles: \n"+v.getProblem());
		// END DEBUGGING

		for (Oracle oracle : oracles) {
			LinkedList<PartialAlignment> segments = oracle
					.findSimiliarSegments(v.getProblem(), M);
			for (PartialAlignment seg : segments) {
				results.add(new RDPSolutionTreeAndNode(v, seg));
			}

		}
		return results.toArray(new RDPSolutionTreeAndNode[0]);
	}

	/**
	 * "The function g_{\vee}: V_{\wedge} \rightarrow V^n_{\vee} generates for
	 * an And-node u \in V_{wedge} a set of n mandatory nodes v_i \in V_{\vee},
	 * which all need to be considered and assembled in order to complete the
	 * current alignment in u."<br />
	 * (From: Protein Threading by Recursive Dynamic Programming. JMB 290,
	 * 757-779)
	 * 
	 * @param u
	 *            the AndNode of the SolutionTree to be worked on
	 * @param t
	 *            the complete SolutionTree
	 * @return subproblems
	 */
	private RDPSolutionTreeOrNode[] gOR(RDPSolutionTreeAndNode u,
			RDPSolutionTree t) {
		// for most oracles there are only 2 subproblems created: To the left
		// and to the right of the new aligned segment.
		LinkedList<RDPSolutionTreeOrNode> results = new LinkedList<RDPSolutionTreeOrNode>();

		LinkedList<RDPProblem> subproblems = u.getPA().getSubProblems();

		// Make new OrNodes
		for (RDPProblem subproblem : subproblems) {
			results.add(new RDPSolutionTreeOrNode(u, subproblem, scoring));
		}

		return results.toArray(new RDPSolutionTreeOrNode[0]);
	}

	/**
	 * filters the given set of AND nodes: removes identical or very similiar
	 * nodes removes alignments contradicting biological and structural
	 * constraints
	 * 
	 * @param uSet
	 */
	private void sf(RDPSolutionTreeAndNode[] uSet) {
		// TODO remove identical nodes

		// TODO remove very similiar nodes

		// TODO remove alignments contradicting biological and structural
		// constraints

	}

	/**
	 * finishes a node
	 * 
	 * @param node
	 *            the node to be finished
	 * @param t
	 *            the tree the node is part of
	 */
	private void finish(RDPSolutionTreeNode node, RDPSolutionTree t) {

		// checkFinal() checks if the node can be finished
		if (node.checkFinal()) {
			if (node instanceof RDPSolutionTreeOrNode) {
				for (RDPSolutionTreeNode child : node.getChilds()) {
					((RDPSolutionTreeOrNode) node)
							.addTAs(((RDPSolutionTreeAndNode) child).getTA());
				}
			} else if (node instanceof RDPSolutionTreeAndNode) {
				if (node.isLeaf()) {
					((RDPSolutionTreeAndNode) node).addTA(new TreeAlignment(
							((RDPSolutionTreeAndNode) node).getPA()));
				} else {
					for (RDPSolutionTreeNode child : node.getChilds()) {
						if (node.ta.isEmpty()) {
							// node.addTAs(child.getTA()); // They are already
							// merged
							// WRONG!
							for (TreeAlignment ta1 : child.getTA()) {
								node.addTA(new TreeAlignment(
										HubeRDP.mergePaT(
												((RDPSolutionTreeAndNode) node)
														.getPA(), ta1)));
							}
						} else {
							LinkedList<TreeAlignment> newTAs = new LinkedList<TreeAlignment>();
							for (TreeAlignment ta1 : child.getTA()) {
								for (TreeAlignment ta2 : node.getTA()) {
									newTAs.add(TreeAlignment.merge(ta1, ta2));
								}
							}
							node.setTA(newTAs);
						}
					}
				}
			}
			node.setFinished(true);

			// finish parent node (if exists)
			if (node.getParent() != null) {
				finish(node.getParent(), t);
			}
		}
	}

	/**
	 * merges a Problem and an Alignment
	 * 
	 * @param prob
	 *            a RDPPloblem
	 * @param ali
	 *            An SequenceAlignment for this Problem
	 * @return a with all the possible PartialAlignments
	 */
	/*
	 * static public TreeAlignment mergePaA(RDPProblem prob, SequenceAlignment
	 * ali) {
	 * 
	 * TreeAlignment result = null;
	 * 
	 * // map is the map of the newly generated alignment int[][] map =
	 * ali.calcMap(); // get first aligned character of the template
	 * 
	 * // calculate paTarStart, paTarEnd, paTemStart, paTemEnd // these are the
	 * first/last aligned characters of the sequences int i = 0; // paTemStart
	 * is the first aligned character of the template sequence in this
	 * subproblem int paTemStart = prob.templateStart; i=0; while (i <
	 * map[0].length && map[0][i] < 0) { i++; } if (i < map[0].length) {
	 * paTemStart += i; }
	 * 
	 * // paTemEnd is the last aligned character of the template sequence in
	 * this subproblem int paTemEnd = prob.templateEnd; i=0; // noch nicht
	 * aligniert while (i < map[0].length && map[0][map[0].length-i-1] < 0) {
	 * i++; } if (i < map[0].length) { paTemEnd -= i; }
	 * 
	 * // paTarStart is the first aligned character of the target sequence in
	 * this subproblem int paTarStart = prob.targetStart; i=0; while (i <
	 * map[1].length && map[1][i] < 0) { i++; } if (i < map[1].length) {
	 * paTarStart += i; }
	 * 
	 * // paTarEnd is the last aligned character of the template sequence in
	 * this subproblem int paTarEnd = prob.targetEnd; i = 0; while (i <
	 * map[1].length && map[1][map[1].length-i-1] < 0) { i++; } if (i <
	 * map[1].length) { paTarEnd -= i; }
	 * 
	 * // calculate new partial alignment (pa) SequenceAlignment pa = null; if
	 * (prob.alignment == null) { // no pa was calculated yet, this means that
	 * this is the first call. pa = ali; } else { // DONE check this code!
	 * ~huberste 2013-02-11 // merge problem.alignment with alignment map =
	 * prob.alignment.calcMap(); char[] temRow = prob.alignment.getRow(0); //
	 * templateRow char[] tarRow = prob.alignment.getRow(1); // targetRow
	 * 
	 * int paTemStartInd = 0; // position in temRow, later first aligned
	 * character int pos = 0; // position in problem.templateSequence while (pos
	 * < prob.templateStart) { if (temRow[paTemStartInd] != '-') { pos++; }
	 * paTemStartInd++; }
	 * 
	 * int paTarStartInd = 0; // position in tarRow, later first aligned
	 * character pos = 0; // position in problem.targetSequence while (pos <
	 * prob.targetStart) { if (tarRow[paTarStartInd] != '-') { pos++; }
	 * paTarStartInd++; }
	 * 
	 * int paTemEndInd = temRow.length - 1; // position in temRow, later last
	 * aligned character pos = prob.templateSequence.length() - 1; // position
	 * in problem.templateSequence while (pos > prob.templateEnd) { if
	 * (temRow[paTemEndInd] != '-') { pos--; } paTemEndInd--; }
	 * 
	 * int paTarEndInd = tarRow.length - 1; // position in tarRow, later last
	 * aligned character pos = prob.targetSequence.length() - 1; // position in
	 * problem.targetSequence while (pos > prob.targetEnd) { if
	 * (tarRow[paTarEndInd] != '-') { pos--; } paTarEndInd--; }
	 * 
	 * // Error handling 143 if (paTemStartInd != paTarStartInd) {
	 * System.err.println
	 * ("Error 143 in Oracle: Alignment lengths don't match."); } // Error
	 * handling 147 if (paTemEndInd != paTarEndInd) {
	 * System.err.println("Error 147 in Oracle: Alignment lengths don't match."
	 * ); }
	 * 
	 * char[][] newRows = new char[2][]; // paTemStartInd = length of old
	 * alignment before new alignment // alignment.length = length of new
	 * alignment // (problem.alignment.length() - (paTemEndInd+1)) = length of
	 * old alignment behind new alignment newRows[0] = new char[paTemStartInd +
	 * ali.length() + (prob.alignment.length() - (paTemEndInd+1))]; newRows[1] =
	 * new char[newRows[0].length]; // make new Alignment: // copy beginning of
	 * old alignment for (pos = 0; pos < paTemStartInd; pos++) {
	 * newRows[0][pos]=temRow[pos]; } for (pos = 0; pos < paTarStartInd; pos++)
	 * { newRows[1][pos]=tarRow[pos]; } // copy new aligned for (pos = 0; pos <
	 * ali.length(); pos++) { newRows[0][paTemStartInd +
	 * pos]=ali.getRow(0)[pos]; } for (pos = 0; pos < ali.length(); pos++) {
	 * newRows[1][paTarStartInd + pos]=ali.getRow(1)[pos]; } // copy end of old
	 * alignment for (pos = 0; pos < prob.alignment.length() - (paTemEndInd +
	 * 1); pos++) { newRows[0][paTemStartInd + ali.length() +
	 * pos]=temRow[(paTemEndInd + 1) + pos]; } for (pos = 0; pos <
	 * prob.alignment.length() - (paTarEndInd + 1); pos++) {
	 * newRows[1][paTarStartInd + ali.length() + pos]=tarRow[(paTarEndInd + 1) +
	 * pos]; }
	 * 
	 * // merge scores. Here simply adding the scores will do fine. double
	 * newScore = prob.alignment.getScore()+ali.getScore();
	 * 
	 * pa = new SequenceAlignment( prob.templateSequence, prob.targetSequence,
	 * newRows, newScore);
	 * 
	 * }
	 * 
	 * result = new TreeAlignment (prob.templateSequence,
	 * prob.templateStructure, prob.targetSequence, prob.targetStructure, pa,
	 * prob.templateStart, prob.templateEnd, prob.targetStart, prob.targetEnd);
	 * 
	 * return result; }
	 */

	/**
	 * merges a Problem and an TreeAlignment
	 * 
	 * @param p1
	 *            a RDP(sub-)Problem
	 * @param p2
	 *            a TreeAlignment for this (sub-)problem
	 * @return a LinkedList with all the possible PartialAlignments
	 */
	static public PartialAlignment mergePaT(RDPProblem p1, RDPProblem p2) {

		// switch them if needed!
		if (p1.alignment.length() < p2.alignment.length()) {
			RDPProblem temp = p1;
			p1 = p2;
			p2 = temp;
		}

		PartialAlignment result = null;

		// int[] alipos = calculateAlignedPositions(ta.alignment);

		char[][] taAliRows = new char[2][];
		taAliRows[0] = p2.alignment.getRow(0);
		taAliRows[1] = p2.alignment.getRow(1);

		char[][] pAliRows = new char[2][];
		pAliRows[0] = p1.alignment.getRow(0);
		pAliRows[1] = p1.alignment.getRow(1);

		int copyTemSeqBefore = p2.templateStart - p1.templateStart + 1;
		int copyTemRowBefore = 0; // length of the part of the p array that
									// needs to be copied before the ta
									// alignment
		int i = 0;
		for (i = 0; i < p1.alignment.length(); i++) {
			if (pAliRows[0][i] != '-') {
				copyTemSeqBefore--;
				if (copyTemSeqBefore == 0) {
					break;
				}
			}
		}
		copyTemRowBefore = i;

		int copyTemSeqBehind = p1.templateEnd - p2.templateEnd + 1;
		int copyTemRowBehind = 0;

		for (i = 0; i < p1.alignment.length(); i++) {
			if (pAliRows[0][p1.alignment.length() - 1 - i] != '-') {
				copyTemSeqBehind--;
				if (copyTemSeqBehind == 0) {
					break;
				}
			}
		}
		copyTemRowBehind = i;

		int copyTarSeqBefore = p2.targetStart - p1.targetStart + 1;
		int copyTarRowBefore = 0;

		for (i = 0; i < p1.alignment.length(); i++) {
			if (pAliRows[1][i] != '-') {
				copyTarSeqBefore--;
				if (copyTarSeqBefore == 0) {
					break;
				}
			}
		}
		copyTarRowBefore = i;

		int copyTarSeqBehind = p1.targetEnd - p2.targetEnd + 1;
		int copyTarRowBehind = 0;

		for (i = 0; i < p1.alignment.length(); i++) {
			if (pAliRows[1][p1.alignment.length() - 1 - i] != '-') {
				copyTarSeqBehind--;
				if (copyTarSeqBehind == 0) {
					break;
				}
			}
		}
		copyTarRowBehind = i;

		int copyBefore = Math.min(copyTemRowBefore, copyTarRowBefore);
		int copyBehind = Math.min(copyTemRowBehind, copyTarRowBehind);

		char[][] newRows = new char[2][];
		newRows[0] = new char[copyBefore + taAliRows[0].length + copyBehind];
		newRows[1] = new char[newRows[0].length];

		// make new Alignment:
		// copy beginning of old alignment
		for (int pos = 0; pos < copyBefore; pos++) {
			newRows[0][pos] = pAliRows[0][pos]; // template row
			newRows[1][pos] = pAliRows[1][pos]; // target row
		}
		// copy ta.alignment
		for (int pos = 0; pos < taAliRows[0].length; pos++) {
			newRows[0][copyBefore + pos] = taAliRows[0][pos]; // template row
			newRows[1][copyBefore + pos] = taAliRows[1][pos]; // target row
		}
		// copy end of old alignment
		for (int pos = 0; pos < copyBehind; pos++) {
			newRows[0][newRows[0].length - pos - 1] = pAliRows[0][pAliRows[0].length
					- pos - 1]; // template row
			newRows[1][newRows[1].length - pos - 1] = pAliRows[1][pAliRows[1].length
					- pos - 1]; // target row
		}

		// TODO merge score
		double newScore = p2.alignment.getScore() + p1.alignment.getScore();

		result = new PartialAlignment(p1.templateSequence,
				p1.templateStructure, p1.targetSequence, p1.targetStructure,
				new SequenceAlignment(p1.templateSequence, p1.targetSequence,
						newRows, newScore), p1.templateStart, p1.templateEnd,
				p1.targetStart, p1.targetEnd, 0, 0, 0, 0);

		return result;
	}

	/**
	 * @return the scoring
	 */
	public Scoring getScoring() {
		return scoring;
	}

	/**
	 * @param scoring
	 *            the scoring to set
	 */
	public void setScoring(Scoring scoring) {
		this.scoring = scoring;
	}

	public static SequenceAlignment hubeRDPAlign(Sequence template,
			Sequence target) {
		// construct rdp tree
		// System.out.print("Constructing RDP Tree structure...");
		RDPProblem root = new RDPProblem(template, null, target, null, null, 0,
				template.length() - 1, 0, target.length() - 1);
		RDPSolutionTree t = new RDPSolutionTree(root);
		// System.out.println(" done!");

		// construct priority queue
		// System.out.print("constructing Priority Queue...");
		RDPPriorityQueue pq = new RDPPriorityQueue(t.getRoot());
		// System.out.println(" done!");

		// construct RDP
		HubeRDP rdp = new HubeRDP();

		// add oracles
		rdp.addOracle(new TinyOracle());
		// rdp.addOracle(new ManualOracle());

		// set scoring
		rdp.setScoring(new SimpleScoring());
		// rdp.setScoring(new ManualScoring());

		// execute rdp algorithm
		// System.out.println("HubeRDP will now be executed!");
		rdp.rdp(t, pq);
		// System.out.println("HubeRDP was successfully executed!");

		// Solution is now in t.getRoot();
		return t.getRoot().getTA().get(0).alignment;
	}

}

/******************************************************************************
 * "A question that sometimes drives me hazy: * Am I or are the others crazy?" *
 * - Albert Einstein (1879 - 1955) *
 ******************************************************************************/
