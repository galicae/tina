package test;

import util.JMolView;

public class JmolTest {

	public static void main(String[] args) {
		String inlinePdb = "ATOM     31  CB  VAL A   4      14.588  -8.758  20.143";
		JMolView jmol = new JMolView(inlinePdb);
	}
}
