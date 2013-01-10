package bioinfo.proteins.fragm3nt;

import java.text.DecimalFormat;
import java.util.LinkedList;
import java.util.Locale;

public class FragmentCluster {
	private LinkedList<ProteinFragment> fragments = new LinkedList<ProteinFragment>();
	private String name;
	private ProteinFragment centroid;
	public double[][] pssm;

	public void calculateCentroid() {
		int actLength = centroid.getAllResidues().length;
		double[][] newCentroid = new double[actLength][3];
		ProteinFragment curFragment;
		double[] curResidue = new double[3];
		while (!fragments.isEmpty()) {
			curFragment = fragments.pop();
			for (int i = 0; i < actLength; i++) {
				curResidue = curFragment.getResidue(i);
				newCentroid[i][0] += curResidue[0];
				newCentroid[i][1] += curResidue[1];
				newCentroid[i][2] += curResidue[2];
			}
		}
		for (int i = 0; i < actLength; i++) {
			newCentroid[i][0] /= actLength * 1.;
			newCentroid[i][1] /= actLength * 1.;
			newCentroid[i][2] /= actLength * 1.;
		}

		centroid = new ProteinFragment(centroid.getID(), newCentroid,
				centroid.getStartIndex(), actLength);
	}

	public ProteinFragment getCentroid() {
		return centroid;
	}

	public void setCentroid(ProteinFragment centroid) {
		this.centroid = centroid;
	}

	public boolean equals(FragmentCluster o) {
		FragmentCluster f = o;
		if (f.getCentroid().getFragmentLength() != centroid.getFragmentLength())
			return false;
		for (int i = 0; i < centroid.getFragmentLength(); i++) {
			for (int j = 0; j < 3; j++) {
				if (!isInEpsilon(f.getCentroid().getAllResidues()[i][j],
						centroid.getAllResidues()[i][j]))
					return false;
			}
		}
		return true;
	}

	public void add(ProteinFragment f) {
		fragments.add(f);
	}

	public void flush() {
		while (!fragments.isEmpty()) {
			// reset everything
			fragments.getFirst().setClusterIndex(-1);
			// then remove it
			fragments.remove();
		}
	}

	private static final double epsilon = 0.0001d;

	private static boolean isInEpsilon(double a, double b) {
		return (a > (b - epsilon)) && (a < (b + epsilon));
	}

	@Override
	public String toString() {
		int i = 0;
		StringBuilder result = new StringBuilder();
		String tempResult = "";

		 String pattern = "####.###";
		 DecimalFormat myFormatter = new DecimalFormat(pattern);

		for (ProteinFragment f : fragments) {
			result.append("REMARK 500 " + f.getSequence() + "\n");
			result.append("MODEL        " + ++i + "\n");
			for (int j = 0; j < f.getAllResidues().length; j++) {
				tempResult = "ATOM  ##### aaaa+rrr c****i   xxxxxxxxyyyyyyyyzzzzzzzzooooootttttt          eehh\n";
				// coordinate strings
				 String xCoord = myFormatter.format(f.getResidue(j)[0]);
//				String xCoord = String.valueOf(f.getResidue(j)[0]);
				while (xCoord.length() < 8)
					xCoord = " " + xCoord;
				 String yCoord = myFormatter.format(f.getResidue(j)[1]);
//				String yCoord = String.valueOf(f.getResidue(j)[1]);
				while (yCoord.length() < 8)
					yCoord = " " + yCoord;
				 String zCoord = myFormatter.format(f.getResidue(j)[2]);
//				String zCoord = String.valueOf(f.getResidue(j)[2]);
				while (zCoord.length() < 8)
					zCoord = " " + zCoord;
				// atom type
				String atomType = " C  ";
				String element = " C";

				tempResult = tempResult.replace("oooooo", "  1.00");
				tempResult = tempResult.replace("tttttt", "  0.00");
				tempResult = tempResult.replace("h", " ");
				tempResult = tempResult.replace("c", "A");
				tempResult = tempResult.replace("+", " ");
				tempResult = tempResult.replace("i", " ");
				tempResult = tempResult.replace("ee", element);
				tempResult = tempResult.replace("****",
						String.format("%4d", (int) j / 4));
				tempResult = tempResult.replace("rrr", "ALA");
				tempResult = tempResult.replace("#####",
						String.format("%5d", j));
				tempResult = tempResult.replace("aaaa", atomType);
				tempResult = tempResult.replace("xxxxxxxx", xCoord);
				tempResult = tempResult.replace("yyyyyyyy", yCoord);
				tempResult = tempResult.replace("zzzzzzzz", zCoord);
				result.append(tempResult);
			}
			tempResult = tempResult
					+ ("TER " + (f.getAllResidues().length - 1) + "\n");
			result.append("ENDMDL" + "\n");
		}
		return result.toString().replace(",", ".");
	}

	public int getSize() {
		return fragments.size();
	}
	
	public LinkedList<ProteinFragment> getFragments() {
		return this.fragments;
	}
	
	public void setPssm(double[][] matrix) {
		pssm = matrix;
	}
}
