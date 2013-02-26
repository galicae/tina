package bioinfo.proteins.corecluster;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

import stat.evaltree.AEvalLeaf;
import stat.evaltree.AEvalNode;
import stat.evaltree.AEvalSupervisor;
import stat.evaltree.IEvalVertix;

public class MultiCurveEvaluationTree extends
		AEvalSupervisor<MultiCurveDataPoint> {

	@Override
	public void generateTree() {
		String[] theta = { "<40", "40-140", ">140" };
		String[] voro = { "important", "unimportant" };
		String[] type = { "HH", "HE", "EE" };
		String[] length = { "short (<7)", "long (>6)" };

		// leafs containing voronoi area classifications
		List<IEvalVertix<MultiCurveDataPoint>>[] voros = new List[18];
		HashMap<String, Integer> voro_class = new HashMap<String, Integer>();
		voro_class.put("important", 0);
		voro_class.put("unimportant", 1);

		for (int i = 0; i != 18; i++) {
			voros[i] = new ArrayList<IEvalVertix<MultiCurveDataPoint>>();
			for (int j = 0; j < 2; j++) {
				voros[i].add(new AEvalLeaf<MultiCurveDataPoint>(voro[j], "") {
					@Override
					public void evalData(MultiCurveDataPoint data) {
						this.addData(data);
					}
				});
			}
		}

		// nodes containing all theta info
		HashMap<String, Integer> theta_class = new HashMap<String, Integer>();
		theta_class.put("<40", 0);
		theta_class.put("40-140", 1);
		theta_class.put(">140", 2);
		List<IEvalVertix<MultiCurveDataPoint>>[] thetas = new List[6];
		for (int i = 0; i != 6; i++) {
			thetas[i] = new ArrayList<IEvalVertix<MultiCurveDataPoint>>();
			for (int j = 0; j != 3; j++) {
				thetas[i].add(new AEvalNode<MultiCurveDataPoint>(theta[j], "",
						voros[i * 3 + j], voro_class) {
					@Override
					public void evalData(MultiCurveDataPoint data) {
						if (data.getVoronoiArea() < 20)
							children.get(classifications.get("unimportant"))
									.evalData(data);
						else
							children.get(classifications.get("important"))
									.evalData(data);
					}

				});
			}
		}

		// nodes containing all type info
		HashMap<String, Integer> type_class = new HashMap<String, Integer>();
		type_class.put("HH", 0);
		type_class.put("HE", 1);
		type_class.put("EE", 2);
		List<IEvalVertix<MultiCurveDataPoint>>[] types = new List[2];

		for (int i = 0; i < 2; i++) {
			types[i] = new ArrayList<IEvalVertix<MultiCurveDataPoint>>();
			for (int j = 0; j != 3; j++) {
				types[i].add(new AEvalNode<MultiCurveDataPoint>(type[j], "",
						thetas[i*3 + j], theta_class) {

					@Override
					public void evalData(MultiCurveDataPoint data) {
						if (data.getTheta() < 40)
							children.get(0).evalData(data);
						else if (data.getTheta() <= 140)
							children.get(1).evalData(data);
						else
							children.get(2).evalData(data);
					}
				});
			}
		}

		// nodes containing all type info
		HashMap<String, Integer> length_class = new HashMap<String, Integer>();
		length_class.put("short", 0);
		length_class.put("long", 1);
		List<IEvalVertix<MultiCurveDataPoint>> lengths = new ArrayList<IEvalVertix<MultiCurveDataPoint>>();

		for (int j = 0; j != 2; j++) {
			lengths.add(new AEvalNode<MultiCurveDataPoint>(length[j], "",
					types[j], type_class) {

				@Override
				public void evalData(MultiCurveDataPoint data) {
					if (data.getType().equals("HH")) {
						children.get(0).evalData(data);
					} else if (data.getType().equals("HE")) {
						children.get(1).evalData(data);
					} else {
						children.get(2).evalData(data);
					}
				}
			});
		}

		// node containing root selection info
		IEvalVertix<MultiCurveDataPoint> root = new AEvalNode<MultiCurveDataPoint>(
				"whole data", "", lengths, length_class) {

			@Override
			public void evalData(MultiCurveDataPoint data) {
				if (data.getMinLength() < 7) {
					children.get(0).evalData(data);
				} else {
					children.get(1).evalData(data);
				}
			}
		};

		this.tree = root;

	}
}

