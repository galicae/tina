package bioinfo.energy.potential.voroEval;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

import bioinfo.proteins.AminoAcidName;

import stat.evaltree.*;

public class VoroEvalTree extends AEvalSupervisor<VoroEvalDataPoint>{

	private double minCon;
	public VoroEvalTree(double minCon){
		this.minCon = minCon;
	}
	
	@Override
	public void generateTree() {
		
		//leafs containing aminoacid classifications
		List<IEvalVertix<VoroEvalDataPoint>>[] leafs = new List[12];
		HashMap<String,Integer> leaf_class = new HashMap<String,Integer>();
		int cnt = 0;
		for(int i = 0; i != 12; i++){
			leafs[i] = new ArrayList<IEvalVertix<VoroEvalDataPoint>>();
			cnt = 0;
			for(AminoAcidName aa : AminoAcidName.values()){
				leafs[i].add(new AEvalLeaf<VoroEvalDataPoint>(aa.getLongName(),""){
					@Override
					public void evalData(VoroEvalDataPoint data) {
						this.addData(data);
					}
				});
				leaf_class.put(aa.getLongName(), cnt);
				cnt++;
			}
		}
		
		//nodes containing all diff info
		String[] diff_desc = {"unterschied<20","20<=unterschied<70","70<unterschied"};
		HashMap<String, Integer> diff_class = new HashMap<String,Integer>();
		diff_class.put("<20",0);
		diff_class.put("20-70",1);
		diff_class.put(">70", 2);
		List<IEvalVertix<VoroEvalDataPoint>>[] diffs = new List[4];
		for(int i = 0; i != 4; i++){
			diffs[i] = new ArrayList<IEvalVertix<VoroEvalDataPoint>>();
			for(int j = 0; j != 3; j++){
				diffs[i].add(new AEvalNode<VoroEvalDataPoint>(diff_desc[j],"",leafs[i*3+j],leaf_class){
								
					@Override
					public void evalData(VoroEvalDataPoint data) {
						children.get(classifications.get(data.getAmino().getLongName())).evalData(data);
					}
					
				});
			}
		}
		
		//nodes containing all dssp inside outside info
		String[] dssp_desc = {"dssp: nicht zugänglich","dssp: zugänglich"};
		HashMap<String, Integer> dssp_class = new HashMap<String,Integer>();
		dssp_class.put("<"+minCon,0);
		dssp_class.put(">="+minCon,1);
		List<IEvalVertix<VoroEvalDataPoint>>[] dssps = new List[2];
		for(int i = 0; i != 2; i++){
			dssps[i] = new ArrayList<IEvalVertix<VoroEvalDataPoint>>();
			for(int j = 0; j != 2; j++){
				dssps[i].add(new AEvalNode<VoroEvalDataPoint>(dssp_desc[j],"",diffs[i*2+j],diff_class) {
					
					@Override
					public void evalData(VoroEvalDataPoint data) {
						if(Math.abs(data.getDsspAcc()-data.getSelfAcc())<20){
							children.get(0).evalData(data);
						}else if(Math.abs(data.getDsspAcc()-data.getSelfAcc())>70){
							children.get(2).evalData(data);
						}else{
							children.get(1).evalData(data);
						}
					}
				});
			}
		}
		
		//nodes containing all self inside outside info
		String[] self_desc = {"self: nicht zugänglich","self: zugänglich"};
		List<IEvalVertix<VoroEvalDataPoint>> selfs;
		HashMap<String, Integer> self_class = new HashMap<String,Integer>();
		diff_class.put("<"+minCon,0);
		diff_class.put(">="+minCon,1);
		selfs = new ArrayList<IEvalVertix<VoroEvalDataPoint>>();
		for(int j = 0; j != 2; j++){
			selfs.add(new AEvalNode<VoroEvalDataPoint>(self_desc[j],"",dssps[j],dssp_class) {
				
				@Override
				public void evalData(VoroEvalDataPoint data) {
					if(data.getDsspAcc()<minCon){
						children.get(0).evalData(data);
					}else{
						children.get(1).evalData(data);
					}
				}
			});
		}
		
		
		//node containing root selection info
		IEvalVertix<VoroEvalDataPoint> root = new AEvalNode<VoroEvalDataPoint>("gesamte Daten","",selfs,self_class){
				
			@Override
			public void evalData(VoroEvalDataPoint data) {
				if(data.getSelfAcc()<minCon){
					children.get(0).evalData(data);
				}else{
					children.get(1).evalData(data);
				}
			}
		};
		
		this.tree = root;
	}
	
}


