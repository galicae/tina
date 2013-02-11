/******************************************************************************
 * hubeRDP is a RDP implementation for the GoBi 2012/13.                      *
 *                                                                            *
 * This file is best read at line width 80 and tab width 4.                   *
 *                                                                   huberste *
 ******************************************************************************/
package huberdp;

import java.util.LinkedList;

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
//					if (u.isLeaf()) { // TODO: u is always leaf?
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
		
		// FIXME
		// BEGIN DEBUGGING
		System.out.println("sending problem to oracles: \n"+v.getProblem());
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
						new RDPSolution(
							((RDPSolutionTreeAndNode) node).getPA().templateSequence,
							((RDPSolutionTreeAndNode) node).getPA().templateStructure,
							((RDPSolutionTreeAndNode) node).getPA().targetSequence,
							((RDPSolutionTreeAndNode) node).getPA().targetStructure,
							((RDPSolutionTreeAndNode) node).getPA().alignment
						)
					);
				} else {
					for (RDPSolutionTreeNode child : node.getChilds()) {
						mergeTAs((RDPSolutionTreeAndNode) node, (RDPSolutionTreeOrNode) child);
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
	
	private PartialAlignment mergeTAs(RDPSolutionTreeAndNode u, RDPSolutionTreeOrNode v) {
		// TODO
		return null;
	}
	
}

/******************************************************************************
 * "A question that sometimes drives me hazy:                                 *
 *  Am I or are the others crazy?"                                            *
 *     - Albert Einstein (1979 - 1966)                                        *
 ******************************************************************************/