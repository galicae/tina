/******************************************************************************
 * huberdp.HubeRDP.java                                                       *
 *                                                                            *
 * This file contains the class HubeRDP which is an RDP implementation for    *
 * the GoBi 2012/13.                                                          *
 *                                                                            *
 * This file is best read at line width 80 and tab width 4.                   *
 *                                                                   huberste *
 ******************************************************************************/
package huberdp;

import java.util.LinkedList;

/**
 * @author huberste
 * @lastchange 2013-02-24
 */
public class HubeRDP {

	private final static int M = 3; // number of different solutions an oracle
									// shall give

	private LinkedList<Oracle> oracles;

	private Scoring scoring;

	/**
	 * constructs an HubeRDP
	 */
	public HubeRDP() {
		oracles = new LinkedList<Oracle>();
	}

	public HubeRDP(Oracle oracle, Scoring scoring) {
		this();
		oracles.add(oracle);
		this.scoring = scoring;
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
	 * @param scoring
	 *            the scoring to set
	 */
	public void setScoring(Scoring scoring) {
		this.scoring = scoring;
	}

	/**
	 * rdp main function. directly transcriped from the RDP paper.
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
			RDPSolutionTreeAndNode[] uSet = gAND(v, M);
			// U <-- sf_{\wedge}(U)
			sf(uSet);
			// if (U = {}) do <SP, TA>^{\wedge} <-- Finish(v,T)
			if (uSet.length == 0) {
				finish(v);
			} else {
				// T <-- insert (T,U)
				v.addChildren(uSet);
				// for each u := <PA, {}>^{\wedge} \in U do
				for (RDPSolutionTreeAndNode u : uSet) {
					// if (Leaf(u)) do <PA,TA>^{\wedge} <-- Finish (u, T)
					// this would always be correct, something's wrong here.
					// Either it's me or the paper. Or my implementation.
					// if (u.isLeaf()) {
					// finish(u,t);
					// } else {
					// V:= {<SP', {}>^{\vee}} <-- g_{\vee}(u, T)
					RDPSolutionTreeOrNode[] vSet = gOR(u);

					if (vSet.length == 0) {
						finish(u);
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
	 * @return (partial) solutions for this subproblem
	 */
	private RDPSolutionTreeAndNode[] gAND(RDPSolutionTreeOrNode v, int m) {

		LinkedList<RDPSolutionTreeAndNode> results = new LinkedList<RDPSolutionTreeAndNode>();

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
	 * @return subproblems
	 */
	private RDPSolutionTreeOrNode[] gOR(RDPSolutionTreeAndNode u) {
		// for most oracles there are only 2 subproblems created: To the left
		// and to the right of the new aligned segment.
		LinkedList<RDPSolutionTreeOrNode> results = new LinkedList<RDPSolutionTreeOrNode>();

		LinkedList<RDPProblem> subproblems = u.getPA().getSubProblems(scoring);

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
	 */
	private void finish(RDPSolutionTreeNode node) {

		// checkFinal() checks if the node can be finished
		if (node.checkFinal()) {
			if (node instanceof RDPSolutionTreeOrNode) {
				for (RDPSolutionTreeNode child : node.getChilds()) {
					((RDPSolutionTreeOrNode) node)
							.addTAs(((RDPSolutionTreeAndNode) child).getTA());
				}
			} else if (node instanceof RDPSolutionTreeAndNode) {
				((RDPSolutionTreeAndNode) node).MergeTAs(scoring);
			}
			node.setFinished(true);

			// finish parent node (if exists)
			if (node.getParent() != null) {
				finish(node.getParent());
			}
		}
	}

}

/******************************************************************************
 * "A question that sometimes drives me hazy:                                 *
 *  Am I or are the others crazy?"                                            *
 * - Albert Einstein (1879 - 1955)                                            *
 ******************************************************************************/
