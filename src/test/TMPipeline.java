package test;

import bioinfo.superpos.TMMain;

public class TMPipeline {
	public static void main(String[] args) throws Exception {
		TMMain main = new TMMain();
		main.calculateTMScore(args[0], args[1], args[2]);
	}
}
