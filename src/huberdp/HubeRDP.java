/******************************************************************************
 * hubeRDP is a RDP implementation for the GoBi 2012/13.                      *
 *                                                                            *
 * This file is best read at line width 80 and tab width 4.                   *
 *                                                                   huberste *
 ******************************************************************************/
package huberdp;

/**
 * @author huberste
 * @lastchange 2013-01-10
 */
public class HubeRDP implements RDP {

	/** 
	 * first rdp must be called with
	 * t = new RDPSolutionTree();
	 * pq = new RDPPriorityQueue(t.getRoot());
	 * rdp (t, pq);
	 */
	@Override
	public void rdp(RDPSolutionTree t, RDPPriorityQueue pq) {
		// TODO Auto-generated method stub
		if (pq.isEmpty() ) {
			return t.getRoot();
		} else {
			// v := <SP, \empty >^{\vee} \leftarrow first(pq)
			RDPSolutionTreeNode v = new RDPSolutionTreeOrNode(pq.getFirst());
			// U := {<PA, \empty >^{\wedge}}
			RDPSolutionTreeNode u = new RDPSolutionTreeAndNode();
			// U <-- sf_{\wedge}(U)
			// TODO
			if (u.isEmpty()) {
				finish(v,t);
			} else {
				
			}
			
		}
	}
	
}

/******************************************************************************
 * "A question that sometimes drives me hazy:                                 *
 *  Am I or are the others crazy?"                                            *
 *     - Albert Einstein (1979 - 1966)                                        *
 ******************************************************************************/