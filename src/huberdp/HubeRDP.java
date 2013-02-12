/******************************************************************************
 * hubeRDP is a RDP implementation for the GoBi 2012/13.                      *
 *                                                                            *
 * This file is best read at line width 80 and tab width 4.                   *
 *                                                                   huberste *
 ******************************************************************************/
package huberdp;

import java.util.LinkedList;

import bioinfo.alignment.SequenceAlignment;

/**
 * @author huberste
 * @lastchange 2013-02-11
 */
public class HubeRDP {

	private final static int M = 3;	// number of different solutions an oracle
									//  shall give
//	private final static int N = 2; // number of subproblems after using oracle
	
	private LinkedList<Oracle> oracles;
	
	/**
	 * constructs an HubeRDP
	 */
	public HubeRDP() {
		oracles = new LinkedList<Oracle>();
	}
	
	/**
	 * adds an oracle to the List of oracles
	 * @param oracle
	 */
	public void addOracle(Oracle oracle) {
		oracles.add(oracle);
	}
	
	/** 
	 * first rdp must be called with
	 * t = new RDPSolutionTree();
	 * pq = new RDPPriorityQueue(t.getRoot());
	 * rdp (t, pq);
	 * Optimal solution is now in t.getRoot();
	 */
	// RDP (T, pq):=
	public void rdp(RDPSolutionTree t, RDPPriorityQueue pq) {
		// if (pq = {} ) do return root
		if (pq.isEmpty() ) {
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
				finish(v,t);
			} else {
				// T <-- insert (T,U)
				v.addChildren(uSet);
				// for each u := <PA, {}>^{\wedge} \in U do
				for (RDPSolutionTreeAndNode u : uSet ) {
					// if (Leaf(u)) do <PA,TA>^{\wedge} <-- Finish (u, T)
					// TODO this would always be correct, something's wrong here
//					if (u.isLeaf()) {
//						finish(u,t);
//					} else {
					// V:= {<SP', {}>^{\vee}} <-- g_{\vee}(u, T)
					RDPSolutionTreeOrNode[] vSet = gOR(u, t);
					
					if (vSet.length == 0) {
						finish(u,t);
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
	 * extensions/completions of the current (partial) alignment coded by 
	 * v."<br />
	 * (From: Protein Threading by Recursive Dynamic Programming. JMB 290,
	 * 757-779)
	 * @param v the OrNode of the SolutionTree to be worked on
	 * @param m maximum number of alternative extensions this function shall
	 * 			create
	 * @param t the complete SolutionTree
	 * @return (partial) solutions for this subproblem
	 */
	private RDPSolutionTreeAndNode[] gAND
			(RDPSolutionTreeOrNode v, int m, RDPSolutionTree t) {
		
		LinkedList<RDPSolutionTreeAndNode> results =
			new LinkedList<RDPSolutionTreeAndNode>();
		
		// BEGIN DEBUGGING
//		System.out.println("sending problem to oracles: \n"+v.getProblem());
		// END DEBUGGING
		
		for (Oracle oracle : oracles) {
		LinkedList<PartialAlignment> segments =
					oracle.findSimiliarSegments(v.getProblem(), M);
			for (PartialAlignment seg: segments) {
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
	 * @param u the AndNode of the SolutionTree to be worked on
	 * @param t the complete SolutionTree
	 * @return subproblems
	 */
	private RDPSolutionTreeOrNode[] gOR
			(RDPSolutionTreeAndNode u, RDPSolutionTree t) {
		// for most oracles there are only 2 subproblems created: To the left
		//  and to the right of the new aligned segment.
		LinkedList<RDPSolutionTreeOrNode> results =
				new LinkedList<RDPSolutionTreeOrNode>();
		
		LinkedList<RDPProblem> subproblems = u.getPA().getSubProblems();
		
		// Make new OrNodes
		for (RDPProblem subproblem : subproblems) {
			results.add( new RDPSolutionTreeOrNode(u, subproblem) );
		}
		
		return results.toArray(new RDPSolutionTreeOrNode[0]);
	}
	
	/**
	 * filters the given set of AND nodes:
	 * removes identical or very similiar nodes
	 * removes alignments contradicting biological and structural constraints
	 * @param uSet
	 */
	private void sf(RDPSolutionTreeAndNode[] uSet) {
		// TODO remove identical nodes
		
		// TODO remove very similiar nodes
		
		// TODO remove alignments contradicting biological and structural
		//	constraints
		
	}
	
	/**
	 * finishes a node
	 * @param node the node to be finished
	 * @param t the tree the node is part of
	 */
	private void finish(RDPSolutionTreeNode node, RDPSolutionTree t) {
		
		// checkFinal() checks if the node can be finished
		if (node.checkFinal() ) {
			if (node instanceof RDPSolutionTreeOrNode) {
				for (RDPSolutionTreeNode child : node.getChilds()) {
					((RDPSolutionTreeOrNode) node).addTAs(((RDPSolutionTreeAndNode) child).getTA());
				}
			} else if (node instanceof RDPSolutionTreeAndNode) {
				if (node.isLeaf()) {
					((RDPSolutionTreeAndNode) node).addTA(
						new TreeAlignment(
							((RDPSolutionTreeAndNode) node).getPA().templateSequence,
							((RDPSolutionTreeAndNode) node).getPA().templateStructure,
							((RDPSolutionTreeAndNode) node).getPA().targetSequence,
							((RDPSolutionTreeAndNode) node).getPA().targetStructure,
							((RDPSolutionTreeAndNode) node).getPA().alignment
						)
					);
				} else {
					for (RDPSolutionTreeNode child : node.getChilds()) {
						for (TreeAlignment solution : child.getTA()) {
							// TODO Look in debug mode what we have right here!
							
							// TODO FIXME
//							LinkedList<PartialAlignment> temp = mergePaA(((RDPSolutionTreeAndNode) node).getPA(), (SequenceAlignment) solution.solution);
//							node.getTA() , solution.solution);
						}
					}
				}
			}
			node.setFinished(true);
			
			// finish parent node (if exists)
			if (node.getParent() != null) {
				finish (node.getParent(), t);
			}
		}
	}
	
	/**
	 * merges a Problem and an Alignment
	 * @param prob a RDPPloblem
	 * @param ali An SequenceAlignment for this Problem
	 * @return a LinkedList with all the possible PartialAlignments
	 */
	static public LinkedList<PartialAlignment> mergePaA(RDPProblem prob, SequenceAlignment ali) {
		
		LinkedList<PartialAlignment> results = new LinkedList<PartialAlignment>();
		
		// map is the map of the newly generated alignment
		int[][] map = ali.calcMap();
		// get first aligned character of the template
		
		// calculate paTarStart, paTarEnd, paTemStart, paTemEnd
		// these are the first/last aligned characters of the sequences
		int i = 0;
		// paTemStart is the first aligned character of the template sequence in this subproblem
		int paTemStart = prob.templateStart;
		i=0;
		while (i < map[0].length && map[0][i] < 0) {
			i++;
		}
		if (i < map[0].length) {
			paTemStart += i;
		}
		
		// paTemEnd is the last aligned character of the template sequence in this subproblem
		int paTemEnd = prob.templateEnd;
		i=0;
		while (i < map[0].length && map[0][map[0].length-i-1] < 0) {
			i++;
		}
		if (i < map[0].length) {
			paTemEnd -= i;
		}
		
		// paTarStart is the first aligned character of the target sequence in this subproblem
		int paTarStart = prob.targetStart;
		i=0;
		while (i < map[1].length && map[1][i] < 0) {
			i++;
		}
		if (i < map[1].length) {
			paTarStart += i;
		}
		
		// paTarEnd is the last aligned character of the template sequence in this subproblem
		int paTarEnd = prob.targetEnd;
		i = 0;
		while (i < map[1].length && map[1][map[1].length-i-1] < 0) {
			i++;
		}
		if (i < map[1].length) {
			paTarEnd -= i;
		}
		
		// calculate new partial alignment (pa)
		SequenceAlignment pa = null;
		if (prob.alignment == null) {
			// no pa was calculated yet, this means that this is the first call.
			pa = ali;
		} else {
			// DONE check this code! ~huberste 2013-02-11
			// merge problem.alignment with alignment
			map = prob.alignment.calcMap();
			char[] temRow = prob.alignment.getRow(0); // templateRow
			char[] tarRow = prob.alignment.getRow(1); // targetRow
			
			int paTemStartInd = 0; // position in temRow, later first aligned character
			int pos = 0; // position in problem.templateSequence
			while (pos < prob.templateStart) {
				if (temRow[paTemStartInd] != '-') {
					pos++;
				}
				paTemStartInd++;
			}
			
			int paTarStartInd = 0; // position in tarRow, later first aligned character
			pos = 0; // position in problem.targetSequence
			while (pos < prob.targetStart) {
				if (tarRow[paTarStartInd] != '-') {
					pos++;
				}
				paTarStartInd++;
			}
			
			int paTemEndInd = temRow.length - 1; // position in temRow, later last aligned character
			pos = prob.templateSequence.length() - 1; // position in problem.templateSequence
			while (pos > prob.templateEnd) {
				if (temRow[paTemEndInd] != '-') {
					pos--;
				}
				paTemEndInd--;
			}
			
			int paTarEndInd = tarRow.length - 1; // position in tarRow, later last aligned character
			pos = prob.targetSequence.length() - 1;  // position in problem.targetSequence
			while (pos > prob.targetEnd) {
				if (tarRow[paTarEndInd] != '-') {
					pos--;
				}
				paTarEndInd--;
			}
			
			// Error handling 143
			if (paTemStartInd != paTarStartInd) {
				System.err.println("Error 143 in Oracle: Alignment lengths don't match.");
			}
			// Error handling 147
			if (paTemEndInd != paTarEndInd) {
				System.err.println("Error 147 in Oracle: Alignment lengths don't match.");
			}
			
			char[][] newRows = new char[2][];
			// paTemStartInd = length of old alignment before new alignment
			// alignment.length = length of new alignment
			// (problem.alignment.length() - (paTemEndInd+1)) = length of old alignment behind new alignment
			newRows[0] = new char[paTemStartInd + ali.length() + (prob.alignment.length() - (paTemEndInd+1))];
			newRows[1] = new char[newRows[0].length];
			// make new Alignment:
			// copy beginning of old alignment
			for (pos = 0; pos < paTemStartInd; pos++) {
				newRows[0][pos]=temRow[pos];
			}
			for (pos = 0; pos < paTarStartInd; pos++) {
				newRows[1][pos]=tarRow[pos];
			}
			// copy new aligned
			for (pos = 0; pos < ali.length(); pos++) {
				newRows[0][paTemStartInd + pos]=ali.getRow(0)[pos];
			}
			for (pos = 0; pos < ali.length(); pos++) {
				newRows[1][paTarStartInd + pos]=ali.getRow(1)[pos];
			}
			// copy end of old alignment
			for (pos = 0; pos < prob.alignment.length() - (paTemEndInd + 1); pos++) {
				newRows[0][paTemStartInd + ali.length() + pos]=temRow[(paTemEndInd + 1) + pos];
			}
			for (pos = 0; pos < prob.alignment.length() - (paTarEndInd + 1); pos++) {
				newRows[1][paTarStartInd + ali.length() + pos]=tarRow[(paTarEndInd + 1) + pos];
			}
			
			// merge scores. Here simply adding the scores will do fine.
			double newScore = prob.alignment.getScore()+ali.getScore();
			
			pa = new SequenceAlignment(
					prob.templateSequence, prob.targetSequence, newRows, newScore);
			
		}
		
		results.add(new PartialAlignment
				(prob.templateSequence, prob.templateStructure,
				prob.targetSequence, prob.targetStructure,
				pa,
				prob.templateStart, prob.templateEnd,
				prob.targetStart, prob.targetEnd,
				paTemStart, paTemEnd, paTarStart, paTarEnd));
		
		return results;
	}
	
}

/******************************************************************************
 * "A question that sometimes drives me hazy:                                 *
 *  Am I or are the others crazy?"                                            *
 *     - Albert Einstein (1979 - 1966)                                        *
 ******************************************************************************/