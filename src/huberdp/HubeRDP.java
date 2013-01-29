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
 * @lastchange 2013-01-26
 */
public class HubeRDP {

	private final static int M = 3;	// number of different solutions an oracle
									//  shall give
//	private final static int N = 2; // number of subproblems after using oracle
	
	private LinkedList<Oracle> oracles;
	
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
					if (u.isLeaf()) { // TODO: u is always leaf?
						finish(u,t);
					} else {
						// V:= {<SP', {}>^{\vee}} <-- g_{\vee}(u, T)
						RDPSolutionTreeOrNode[] vSet = gOR(u, t);
						// T <-- insert (T, V)
						u.addChildren(vSet);
						// pq <-- insert(pq, V)
						pq.add(vSet);
					}
					
				}
				
			}
//			return rdp(t, pq);
		}
	}
	
	/**
	 * Function that calls the oracles
	 * @param v the OrNode of the SolutionTree to be worked on
	 * @param m maximum number of "hits" an oracle shall create
	 * @param t the complete SolutionTree
	 * @return (partial) solutions for this subproblem
	 */
	private RDPSolutionTreeAndNode[] gAND
			(RDPSolutionTreeOrNode v, int m, RDPSolutionTree t) {
		
		LinkedList<RDPSolutionTreeAndNode> results =
			new LinkedList<RDPSolutionTreeAndNode>();
		
		for (Oracle oracle : oracles) {
			LinkedList<RDPProblem> temp =
					oracle.findSimiliarSegments(v.getProblem(), M);
			for (RDPProblem node: temp) {
				results.add(node);
			}
			
		}
		return results.toArray(new RDPSolutionTreeAndNode[0]);
	}
	
	/**
	 * 
	 * @param u the AndNode of the SolutionTree to be worked on
	 * @param t the complete SolutionTree
	 * @return subproblems
	 */
	private RDPSolutionTreeOrNode[] gOR
			(RDPSolutionTreeAndNode u, RDPSolutionTree t) {
		RDPSolutionTreeOrNode[] results = new RDPSolutionTreeOrNode[2];
		// TODO
		u.
		
		return results;
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
		
		// TODO remove alignments contradicting biological and structural constraints
		
	}
	
	
	private void finish(RDPSolutionTreeNode u, RDPSolutionTree t) {
		// TODO
	}
	
}

/******************************************************************************
 * "A question that sometimes drives me hazy:                                 *
 *  Am I or are the others crazy?"                                            *
 *     - Albert Einstein (1979 - 1966)                                        *
 ******************************************************************************/