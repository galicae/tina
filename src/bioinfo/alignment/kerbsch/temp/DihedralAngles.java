package bioinfo.alignment.kerbsch.temp;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;

import bioinfo.alignment.gotoh.Gotoh;

public class DihedralAngles {
	private static final ArrayList<Character> aminoacids = new ArrayList<Character>() {
		{
			add('A');
			add('C');
			add('D');
			add('E');
			add('F');
			add('G');
			add('H');
			add('I');
			add('K');
			add('L');
			add('M');
			add('N');
			add('P');
			add('Q');
			add('R');
			add('S');
			add('T');
			add('V');
			add('W');
			add('Y');
		}
	};

	public static void calc(String folderpath) {
		BufferedReader in;
		BufferedWriter out;

		double[][][][] minusminus = new double[26][26][26][3];
		double[][][][] minusplus = new double[26][26][26][3];
		double[][][][] plusminus = new double[26][26][26][3];
		double[][][][] plusplus = new double[26][26][26][3];

		char precursor;
		char residue;
		char follower;
		double psi = 0.0;
		double phi = 0.0;

		File dir = new File(folderpath);
		File[] files = dir.listFiles();
		int index = 0;
		for (File f : files) {
			System.out.println(++index);
			try {
				in = new BufferedReader(new FileReader(f));
				String line;
				while ((line = in.readLine()) != null) {
					if (line.contains("NUMBER OF RESIDUES")) {
						if (Integer.parseInt(line.split("\\s+")[1]) < 3) {
							break;
						}
					}

					// read angles for every residue except of first and last
					// residue
					// and store them in HashMap
					else if (line.trim().startsWith("#")) {
						System.out.println(f.getName());

						// init steps
						line = in.readLine();
						if (!line.matches("(.)*!(.)*")
								&& line.charAt(13) != 'X') { // chain border
							precursor = line.charAt(13);

							// lower-case letters in DSSP are Cysteines
							if (Character.isLowerCase(precursor)) {
								precursor = 'C';
							}

						} else {
							precursor = 0;
						}
						line = in.readLine();
						if (!line.matches("(.)*!(.)*")
								&& line.charAt(13) != 'X') { // chain border
							residue = line.charAt(13);

							// lower-case letters in DSSP are Cysteines
							if (Character.isLowerCase(residue)) {
								residue = 'C';
							}

							phi = Double.parseDouble(line.substring(103, 109));
							psi = Double.parseDouble(line.substring(109, 115));
						} else {
							residue = 0;
						}
						line = in.readLine();
						if (!line.matches("(.)*!(.)*")
								&& line.charAt(13) != 'X') { // chain border
							follower = line.charAt(13);

							// lower-case letters in DSSP are Cysteines
							if (Character.isLowerCase(follower)) {
								follower = 'C';
							}

						} else {
							follower = 0;
						}
						if (precursor != 0 && residue != 0 && follower != 0) {
							if (phi < 0 && psi < 0) {
								minusminus[residue - 65][precursor - 65][follower - 65][0] += phi;
								minusminus[residue - 65][precursor - 65][follower - 65][1] += psi;
								minusminus[residue - 65][precursor - 65][follower - 65][2] += 1;
							} else if (phi < 0 && psi >= 0) {
								minusplus[residue - 65][precursor - 65][follower - 65][0] += phi;
								minusplus[residue - 65][precursor - 65][follower - 65][1] += psi;
								minusplus[residue - 65][precursor - 65][follower - 65][2] += 1;
							} else if (phi >= 0 && psi < 0) {
								plusminus[residue - 65][precursor - 65][follower - 65][0] += phi;
								plusminus[residue - 65][precursor - 65][follower - 65][1] += psi;
								plusminus[residue - 65][precursor - 65][follower - 65][2] += 1;
							} else if (phi >= 0 && psi >= 0) {
								plusplus[residue - 65][precursor - 65][follower - 65][0] += phi;
								plusplus[residue - 65][precursor - 65][follower - 65][1] += psi;
								plusplus[residue - 65][precursor - 65][follower - 65][2] += 1;
							}
							precursor = residue;
							residue = follower;
							phi = Double.parseDouble(line.substring(103, 109));
							psi = Double.parseDouble(line.substring(109, 115));
						} else {
							precursor = residue;
							residue = follower;
						}

						// read the rest of the DSSP
						while ((line = in.readLine()) != null) {
							if (!line.matches("(.)*!(.)*")
									&& line.charAt(13) != 'X') { // chain border
								follower = line.charAt(13);

								// lower-case letters in DSSP are Cysteines
								if (Character.isLowerCase(follower)) {
									follower = 'C';
								}

								if (precursor != 0 && residue != 0) {
									if (phi < 0 && psi < 0) {
										minusminus[residue - 65][precursor - 65][follower - 65][0] += phi;
										minusminus[residue - 65][precursor - 65][follower - 65][1] += psi;
										minusminus[residue - 65][precursor - 65][follower - 65][2] += 1;
									} else if (phi < 0 && psi >= 0) {
										minusplus[residue - 65][precursor - 65][follower - 65][0] += phi;
										minusplus[residue - 65][precursor - 65][follower - 65][1] += psi;
										minusplus[residue - 65][precursor - 65][follower - 65][2] += 1;
									} else if (phi >= 0 && psi < 0) {
										plusminus[residue - 65][precursor - 65][follower - 65][0] += phi;
										plusminus[residue - 65][precursor - 65][follower - 65][1] += psi;
										plusminus[residue - 65][precursor - 65][follower - 65][2] += 1;
									} else if (phi >= 0 && psi >= 0) {
										plusplus[residue - 65][precursor - 65][follower - 65][0] += phi;
										plusplus[residue - 65][precursor - 65][follower - 65][1] += psi;
										plusplus[residue - 65][precursor - 65][follower - 65][2] += 1;
									}
									precursor = residue;
									residue = follower;
									phi = Double.parseDouble(line.substring(
											103, 109));
									psi = Double.parseDouble(line.substring(
											109, 115));
								} else {
									precursor = residue;
									residue = follower;
								}
							} else {
								precursor = residue;
								residue = follower;
							}
						}
					}
				}
				in.close();
			} catch (IOException e) {
				System.out.println("No DSSP input!");
			}
		}

		// write out results
		for (char amino : aminoacids) {
			try {
				out = new BufferedWriter(new FileWriter("angles/" + amino
						+ ".angles"));
				for (char pre : aminoacids) {
					for (char fol : aminoacids) {
						out.append(pre + "" + fol + "\t");
						
						//minus minus
						if(minusminus[amino - 65][pre - 65][fol - 65][2] > 0){
						out.append(minusminus[amino - 65][pre - 65][fol - 65][0]
								/ minusminus[amino - 65][pre - 65][fol - 65][2]
								+ "\t");
						out.append(minusminus[amino - 65][pre - 65][fol - 65][1]
								/ minusminus[amino - 65][pre - 65][fol - 65][2]
								+ "\t");
						} else{
							out.append(-999.0 + "\t"); //profile not exists
							out.append(-999.0 + "\t");	
						}
						
						//minus plus
						if(minusplus[amino - 65][pre - 65][fol - 65][2] > 0){
						out.append(minusplus[amino - 65][pre - 65][fol - 65][0]
								/ minusplus[amino - 65][pre - 65][fol - 65][2]
								+ "\t");
						out.append(minusplus[amino - 65][pre - 65][fol - 65][1]
								/ minusplus[amino - 65][pre - 65][fol - 65][2]
								+ "\t");
						} else{
							out.append(-999.0 + "\t"); //profile not exists
							out.append(999.0 + "\t");	
						}
						
						
						//plus minus
						if(plusminus[amino - 65][pre - 65][fol - 65][2] > 0){
						out.append(plusminus[amino - 65][pre - 65][fol - 65][0]
								/ plusminus[amino - 65][pre - 65][fol - 65][2]
								+ "\t");
						out.append(plusminus[amino - 65][pre - 65][fol - 65][1]
								/ plusminus[amino - 65][pre - 65][fol - 65][2]
								+ "\t");
						} else{
							out.append(999.0 + "\t"); //profile not exists
							out.append(-999.0 + "\t");	
						}
						
						//plus plus
						if(plusplus[amino - 65][pre - 65][fol - 65][2] > 0){
						out.append(plusplus[amino - 65][pre - 65][fol - 65][0]
								/ plusplus[amino - 65][pre - 65][fol - 65][2]
								+ "\t");
						out.append(plusplus[amino - 65][pre - 65][fol - 65][1]
								/ plusplus[amino - 65][pre - 65][fol - 65][2]
								+ "\n");
						}else{
							out.append(999.0 + "\t"); //profile not exists
							out.append(999.0 + "\t");	
						}
					}
				}
				out.close();
			} catch (IOException e) {
				System.out.println("cannot write angles " + amino);
			}
		}
	}

	public static int[] readStances(String filepath) {
		BufferedReader in;
		String line;
		int[] stances = null;
		double phi;
		double psi;
		int index = -1;

		try {
			in = new BufferedReader(new FileReader(filepath));
			while ((line = in.readLine()) != null) {
				if (line.contains("NUMBER OF RESIDUES")) {
					stances = new int[Integer.parseInt(line.split("\\s+")[1])];
				} else if (line.trim().startsWith("#")) {
					while ((line = in.readLine()) != null) {
						if (!line.matches("(.)*!(.)*") && line.charAt(13) != 'X') { // chain border
							phi = Double.parseDouble(line.substring(103, 109));
							psi = Double.parseDouble(line.substring(109, 115));
	
							if (phi < 0 && psi < 0) {
								stances[++index] = 0;
							} else if (phi < 0 && psi >= 0) {
								stances[++index] = 1;
							} else if (phi >= 0 && psi < 0) {
								stances[++index] = 2;
							} else if (phi >= 0 && psi >= 0) {
								stances[++index] = 3;
							} else {
								stances[++index] = -1; // just in case there is a faulty DSSP
							}
						}
					}
				}
			}
			in.close();
		} catch (IOException e) {
			System.out.println("cannot read Stances " + filepath);
		}
		return stances;
	}

	public static int[][][][][] readAngles(String folderpath) {
		int[][][][][] matrix = new int[26][26][26][2][4];
		File dir = new File(folderpath);
		File[] files = dir.listFiles();
		BufferedReader in;
		String line;
		String[] temp;

		int mmphi;
		int mmpsi;
		int mpphi;
		int mppsi;
		int pmphi;
		int pmpsi;
		int ppphi;
		int pppsi;

		for (File f : files) {
			try {
				in = new BufferedReader(new FileReader(f));
				while ((line = in.readLine()) != null) {
					temp = line.split("\\s+");
					mmphi = (int) (Double.parseDouble(temp[1]) * Gotoh.FACTOR);
					mmpsi = (int) (Double.parseDouble(temp[2]) * Gotoh.FACTOR);
					mpphi = (int) (Double.parseDouble(temp[3]) * Gotoh.FACTOR);
					mppsi = (int) (Double.parseDouble(temp[4]) * Gotoh.FACTOR);
					pmphi = (int) (Double.parseDouble(temp[5]) * Gotoh.FACTOR);
					pmpsi = (int) (Double.parseDouble(temp[6]) * Gotoh.FACTOR);
					ppphi = (int) (Double.parseDouble(temp[7]) * Gotoh.FACTOR);
					pppsi = (int) (Double.parseDouble(temp[8]) * Gotoh.FACTOR);

					matrix[f.getName().charAt(0) - 65][temp[0].charAt(0) - 65][temp[0]
							.charAt(1) - 65][0][0] = mmphi;
					matrix[f.getName().charAt(0) - 65][temp[0].charAt(0) - 65][temp[0]
							.charAt(1) - 65][1][0] = mmpsi;
					matrix[f.getName().charAt(0) - 65][temp[0].charAt(0) - 65][temp[0]
							.charAt(1) - 65][0][1] = mpphi;
					matrix[f.getName().charAt(0) - 65][temp[0].charAt(0) - 65][temp[0]
							.charAt(1) - 65][1][1] = mppsi;
					matrix[f.getName().charAt(0) - 65][temp[0].charAt(0) - 65][temp[0]
							.charAt(1) - 65][0][2] = pmphi;
					matrix[f.getName().charAt(0) - 65][temp[0].charAt(0) - 65][temp[0]
							.charAt(1) - 65][1][2] = pmpsi;
					matrix[f.getName().charAt(0) - 65][temp[0].charAt(0) - 65][temp[0]
							.charAt(1) - 65][0][3] = ppphi;
					matrix[f.getName().charAt(0) - 65][temp[0].charAt(0) - 65][temp[0]
							.charAt(1) - 65][1][3] = pppsi;
				}
				in.close();
			} catch (IOException e) {
				System.out.println("cannot read anglefile " + f);
			}
		}
		return matrix;
	}

	public static void main(String[] args) {
		DihedralAngles.calc("DSSP");
	}

}
