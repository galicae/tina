package bioinfo.proteins;

public enum SecStructThree {

	H('H'), E('E'), C('C');

//	static SecStructThree[] map = new SecStructThree['Z'-'A'];
	private final char charRepres;

	private SecStructThree(char charRepres) {
		if (charRepres == 'L')
			charRepres = 'H';
		this.charRepres = charRepres;
	}

	public char getCharRepres() {
		return this.charRepres;
	}

	public static SecStructThree defSecStructThree(char charRepres) {
		for (SecStructThree ssE : SecStructThree.values()) {
			if (ssE.getCharRepres() == charRepres) {
				return ssE;
			}
		}
		return SecStructThree.C;
	}

}
