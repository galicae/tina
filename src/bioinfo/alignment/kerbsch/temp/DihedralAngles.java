package bioinfo.alignment.kerbsch.temp;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;

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
	private static final HashMap<Character,Integer> aaToIndex = new HashMap<Character,Integer>() {
		{
			put('A',1);
			put('C',2);
			put('D',3);
			put('E',4);
			put('F',5);
			put('G',6);
			put('H',7);
			put('I',8);
			put('K',9);
			put('L',10);
			put('M',11);
			put('N',12);
			put('P',13);
			put('Q',14);
			put('R',15);
			put('S',16);
			put('T',17);
			put('V',18);
			put('W',19);
			put('Y',20);
		}
	};

	public static void calc(String folderpath) {
		BufferedReader in;
		BufferedWriter out;
		HashMap<String,Integer> profilesToIndex = new HashMap<String,Integer>();
		int profindex = 1;
		for(char c : aminoacids){
			for(char c2: aminoacids){
				profilesToIndex.put(c+c2+"",profindex);
				profindex++;
			}
		}
		

		//first dimension: amino acid in the middle
		//second dimension: amino acid precursor
		//third dimension: amino acid follower
		//fourth dimension: 0 -> phi, 1 -> psi, 2 -> count, 3 -> max phi, 4 -> max psi, 5 -> min phi, 6 -> min psi
		double[][][][] minusminus = new double[26][26][26][7];
		double[][][][] minusplus = new double[26][26][26][7];
		double[][][][] plusminus = new double[26][26][26][7];
		double[][][][] plusplus = new double[26][26][26][7];
		
		for (int i = 0; i < 26; i++) {
			for (int j = 0; j < 26; j++) {
				for (int j2 = 0; j2 < 26; j2++) {
					minusminus[i][j][j2][3] = Double.NEGATIVE_INFINITY;
					minusminus[i][j][j2][4] = Double.NEGATIVE_INFINITY;
					minusminus[i][j][j2][5] = Double.POSITIVE_INFINITY;
					minusminus[i][j][j2][6] = Double.POSITIVE_INFINITY;
					
					minusplus[i][j][j2][3] = Double.NEGATIVE_INFINITY;
					minusplus[i][j][j2][4] = Double.NEGATIVE_INFINITY;
					minusplus[i][j][j2][5] = Double.POSITIVE_INFINITY;
					minusplus[i][j][j2][6] = Double.POSITIVE_INFINITY;
					
					plusminus[i][j][j2][3] = Double.NEGATIVE_INFINITY;
					plusminus[i][j][j2][4] = Double.NEGATIVE_INFINITY;
					plusminus[i][j][j2][5] = Double.POSITIVE_INFINITY;
					plusminus[i][j][j2][6] = Double.POSITIVE_INFINITY;
					
					plusplus[i][j][j2][3] = Double.NEGATIVE_INFINITY;
					plusplus[i][j][j2][4] = Double.NEGATIVE_INFINITY;
					plusplus[i][j][j2][5] = Double.POSITIVE_INFINITY;
					plusplus[i][j][j2][6] = Double.POSITIVE_INFINITY;
				}
			}
		}

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
							if(precursor != 0){
								residue = line.charAt(13);
	
								// lower-case letters in DSSP are Cysteines
								if (Character.isLowerCase(residue)) {
									residue = 'C';
								}
								phi = Double.parseDouble(line.substring(103, 109));
								psi = Double.parseDouble(line.substring(109, 115));
							} else{
								residue = 0;
							}
						} else {
							residue = 0;
						}
						line = in.readLine();
						if (!line.matches("(.)*!(.)*")
								&& line.charAt(13) != 'X') { // chain border
							if(residue != 0){
								follower = line.charAt(13);
	
								// lower-case letters in DSSP are Cysteines
								if (Character.isLowerCase(follower)) {
									follower = 'C';
								}
							} else{
								follower = 0;
							}
						} else {
							follower = 0;
						}
						if (precursor != 0 && residue != 0 && follower != 0) {
							if (phi < 0 && psi < 0) {
								//stats | max,min
								if(minusminus[residue - 65][precursor - 65][follower - 65][3] < phi){
									minusminus[residue - 65][precursor - 65][follower - 65][3] = phi;
								}
								if(minusminus[residue - 65][precursor - 65][follower - 65][4] < psi){
									minusminus[residue - 65][precursor - 65][follower - 65][4] = psi;
								}
								if(minusminus[residue - 65][precursor - 65][follower - 65][5] > phi){
									minusminus[residue - 65][precursor - 65][follower - 65][5] = phi;
								}
								if(minusminus[residue - 65][precursor - 65][follower - 65][6] > psi){
									minusminus[residue - 65][precursor - 65][follower - 65][6] = psi;
								}
								minusminus[residue - 65][precursor - 65][follower - 65][0] += phi;
								minusminus[residue - 65][precursor - 65][follower - 65][1] += psi;
								minusminus[residue - 65][precursor - 65][follower - 65][2] += 1;
							} else if (phi < 0 && psi >= 0) {
								//stats | max,min
								if(minusplus[residue - 65][precursor - 65][follower - 65][3] < phi){
									minusplus[residue - 65][precursor - 65][follower - 65][3] = phi;
								}
								if(minusplus[residue - 65][precursor - 65][follower - 65][4] < psi){
									minusplus[residue - 65][precursor - 65][follower - 65][4] = psi;
								}
								if(minusplus[residue - 65][precursor - 65][follower - 65][5] > phi){
									minusplus[residue - 65][precursor - 65][follower - 65][5] = phi;
								}
								if(minusplus[residue - 65][precursor - 65][follower - 65][6] > psi){
									minusplus[residue - 65][precursor - 65][follower - 65][6] = psi;
								}
								minusplus[residue - 65][precursor - 65][follower - 65][0] += phi;
								minusplus[residue - 65][precursor - 65][follower - 65][1] += psi;
								minusplus[residue - 65][precursor - 65][follower - 65][2] += 1;
							} else if (phi >= 0 && psi < 0) {
								//stats | max,min
								if(plusminus[residue - 65][precursor - 65][follower - 65][3] < phi){
									plusminus[residue - 65][precursor - 65][follower - 65][3] = phi;
								}
								if(plusminus[residue - 65][precursor - 65][follower - 65][4] < psi){
									plusminus[residue - 65][precursor - 65][follower - 65][4] = psi;
								}
								if(plusminus[residue - 65][precursor - 65][follower - 65][5] > phi){
									plusminus[residue - 65][precursor - 65][follower - 65][5] = phi;
								}
								if(plusminus[residue - 65][precursor - 65][follower - 65][6] > psi){
									plusminus[residue - 65][precursor - 65][follower - 65][6] = psi;
								}
								plusminus[residue - 65][precursor - 65][follower - 65][0] += phi;
								plusminus[residue - 65][precursor - 65][follower - 65][1] += psi;
								plusminus[residue - 65][precursor - 65][follower - 65][2] += 1;
							} else if (phi >= 0 && psi >= 0) {
								//stats | max,min
								if(plusplus[residue - 65][precursor - 65][follower - 65][3] < phi){
									plusplus[residue - 65][precursor - 65][follower - 65][3] = phi;
								}
								if(plusplus[residue - 65][precursor - 65][follower - 65][4] < psi){
									plusplus[residue - 65][precursor - 65][follower - 65][4] = psi;
								}
								if(plusplus[residue - 65][precursor - 65][follower - 65][5] > phi){
									plusplus[residue - 65][precursor - 65][follower - 65][5] = phi;
								}
								if(plusplus[residue - 65][precursor - 65][follower - 65][6] > psi){
									plusplus[residue - 65][precursor - 65][follower - 65][6] = psi;
								}
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
							phi = Double.parseDouble(line.substring(103, 109));
							psi = Double.parseDouble(line.substring(109, 115));
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
										//stats | max,min
										if(minusminus[residue - 65][precursor - 65][follower - 65][3] < phi){
											minusminus[residue - 65][precursor - 65][follower - 65][3] = phi;
										}
										if(minusminus[residue - 65][precursor - 65][follower - 65][4] < psi){
											minusminus[residue - 65][precursor - 65][follower - 65][4] = psi;
										}
										if(minusminus[residue - 65][precursor - 65][follower - 65][5] > phi){
											minusminus[residue - 65][precursor - 65][follower - 65][5] = phi;
										}
										if(minusminus[residue - 65][precursor - 65][follower - 65][6] > psi){
											minusminus[residue - 65][precursor - 65][follower - 65][6] = psi;
										}
										minusminus[residue - 65][precursor - 65][follower - 65][0] += phi;
										minusminus[residue - 65][precursor - 65][follower - 65][1] += psi;
										minusminus[residue - 65][precursor - 65][follower - 65][2] += 1;
									} else if (phi < 0 && psi >= 0) {
										//stats | max,min
										if(minusplus[residue - 65][precursor - 65][follower - 65][3] < phi){
											minusplus[residue - 65][precursor - 65][follower - 65][3] = phi;
										}
										if(minusplus[residue - 65][precursor - 65][follower - 65][4] < psi){
											minusplus[residue - 65][precursor - 65][follower - 65][4] = psi;
										}
										if(minusplus[residue - 65][precursor - 65][follower - 65][5] > phi){
											minusplus[residue - 65][precursor - 65][follower - 65][5] = phi;
										}
										if(minusplus[residue - 65][precursor - 65][follower - 65][6] > psi){
											minusplus[residue - 65][precursor - 65][follower - 65][6] = psi;
										}
										minusplus[residue - 65][precursor - 65][follower - 65][0] += phi;
										minusplus[residue - 65][precursor - 65][follower - 65][1] += psi;
										minusplus[residue - 65][precursor - 65][follower - 65][2] += 1;
									} else if (phi >= 0 && psi < 0) {
										//stats | max,min
										if(plusminus[residue - 65][precursor - 65][follower - 65][3] < phi){
											plusminus[residue - 65][precursor - 65][follower - 65][3] = phi;
										}
										if(plusminus[residue - 65][precursor - 65][follower - 65][4] < psi){
											plusminus[residue - 65][precursor - 65][follower - 65][4] = psi;
										}
										if(plusminus[residue - 65][precursor - 65][follower - 65][5] > phi){
											plusminus[residue - 65][precursor - 65][follower - 65][5] = phi;
										}
										if(plusminus[residue - 65][precursor - 65][follower - 65][6] > psi){
											plusminus[residue - 65][precursor - 65][follower - 65][6] = psi;
										}
										plusminus[residue - 65][precursor - 65][follower - 65][0] += phi;
										plusminus[residue - 65][precursor - 65][follower - 65][1] += psi;
										plusminus[residue - 65][precursor - 65][follower - 65][2] += 1;
									} else if (phi >= 0 && psi >= 0) {
										//stats | max,min
										if(plusplus[residue - 65][precursor - 65][follower - 65][3] < phi){
											plusplus[residue - 65][precursor - 65][follower - 65][3] = phi;
										}
										if(plusplus[residue - 65][precursor - 65][follower - 65][4] < psi){
											plusplus[residue - 65][precursor - 65][follower - 65][4] = psi;
										}
										if(plusplus[residue - 65][precursor - 65][follower - 65][5] > phi){
											plusplus[residue - 65][precursor - 65][follower - 65][5] = phi;
										}
										if(plusplus[residue - 65][precursor - 65][follower - 65][6] > psi){
											plusplus[residue - 65][precursor - 65][follower - 65][6] = psi;
										}
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
									phi = Double.parseDouble(line.substring(103, 109));
									psi = Double.parseDouble(line.substring(109, 115));
								}
							} else {
								follower = 0;
								precursor = residue;
								residue = follower;
								phi = Double.parseDouble(line.substring(103, 109));
								psi = Double.parseDouble(line.substring(109, 115));
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
		try {
			out = new BufferedWriter(new FileWriter("angles/summary.angles"));
			for (char amino : aminoacids) {
				for (char pre : aminoacids) {
					for (char fol : aminoacids) {
//						out.append(pre + "" + fol + "\t");
						
						out.append(aaToIndex.get(amino)+ "\t" + profilesToIndex.get(pre+fol+"") + "\t");
						
						//minus minus
						if(minusminus[amino - 65][pre - 65][fol - 65][2] > 0){
						out.append(Math.round(minusminus[amino - 65][pre - 65][fol - 65][0]
								/ minusminus[amino - 65][pre - 65][fol - 65][2]*100.0)/100.0
								+ "\t");
						out.append(Math.round(minusminus[amino - 65][pre - 65][fol - 65][1]
								/ minusminus[amino - 65][pre - 65][fol - 65][2]*100.0)/100.0
								+ "\t");
						//counts
						out.append(minusminus[amino - 65][pre - 65][fol - 65][2] + "\t");
						
						//maxmin phi
						out.append(minusminus[amino - 65][pre - 65][fol - 65][3] + "\t");
						out.append(minusminus[amino - 65][pre - 65][fol - 65][5] + "\t");
						
						//maxmin psi
						out.append(minusminus[amino - 65][pre - 65][fol - 65][4] + "\t");
						out.append(minusminus[amino - 65][pre - 65][fol - 65][6] + "\t");
						} else{
							out.append(-999.0 + "\t"); //profile not exists
							out.append(-999.0 + "\t" + 0.0 + "\t");	
							out.append(360.0 + "\t" + -360.0 + "\t" );
							out.append(360.0 + "\t" + -360.0 + "\t" );
						}
						
						//minus plus
						if(minusplus[amino - 65][pre - 65][fol - 65][2] > 0){
						out.append(Math.round(minusplus[amino - 65][pre - 65][fol - 65][0]
								/ minusplus[amino - 65][pre - 65][fol - 65][2]*100.0)/100.0
								+ "\t");
						out.append(Math.round(minusplus[amino - 65][pre - 65][fol - 65][1]
								/ minusplus[amino - 65][pre - 65][fol - 65][2]*100.0)/100.0
								+ "\t");
						//counts
						out.append(minusplus[amino - 65][pre - 65][fol - 65][2] + "\t");
						
						//maxmin phi
						out.append(minusplus[amino - 65][pre - 65][fol - 65][3] + "\t");
						out.append(minusplus[amino - 65][pre - 65][fol - 65][5] + "\t");
						
						//maxmin psi
						out.append(minusplus[amino - 65][pre - 65][fol - 65][4] + "\t");
						out.append(minusplus[amino - 65][pre - 65][fol - 65][6] + "\t");
						} else{
							out.append(-999.0 + "\t"); //profile not exists
							out.append(999.0 + "\t" + 0.0 + "\t");	
							out.append(360.0 + "\t" + -360.0 + "\t" );
							out.append(360.0 + "\t" + -360.0 + "\t" );
						}
						
						
						//plus minus
						if(plusminus[amino - 65][pre - 65][fol - 65][2] > 0){
						out.append(Math.round(plusminus[amino - 65][pre - 65][fol - 65][0]
								/ plusminus[amino - 65][pre - 65][fol - 65][2]*100.0)/100.0
								+ "\t");
						out.append(Math.round(plusminus[amino - 65][pre - 65][fol - 65][1]
								/ plusminus[amino - 65][pre - 65][fol - 65][2]*100.0)/100.0
								+ "\t");
						//counts
						out.append(plusminus[amino - 65][pre - 65][fol - 65][2] + "\t");
						
						//maxmin phi
						out.append(plusminus[amino - 65][pre - 65][fol - 65][3] + "\t");
						out.append(plusminus[amino - 65][pre - 65][fol - 65][5] + "\t");
						
						//maxmin psi
						out.append(plusminus[amino - 65][pre - 65][fol - 65][4] + "\t");
						out.append(plusminus[amino - 65][pre - 65][fol - 65][6] + "\t");
						} else{
							out.append(999.0 + "\t"); //profile not exists
							out.append(-999.0 + "\t" + 0.0 + "\t");	
							out.append(360.0 + "\t" + -360.0 + "\t" );
							out.append(360.0 + "\t" + -360.0 + "\t" );
						}
						
						//plus plus
						if(plusplus[amino - 65][pre - 65][fol - 65][2] > 0){
						out.append(Math.round(plusplus[amino - 65][pre - 65][fol - 65][0]
								/ plusplus[amino - 65][pre - 65][fol - 65][2]*100.0)/100.0
								+ "\t");
						out.append(Math.round(plusplus[amino - 65][pre - 65][fol - 65][1]
								/ plusplus[amino - 65][pre - 65][fol - 65][2]*100.0)/100.0
								+ "\t");
						//counts
						out.append(plusplus[amino - 65][pre - 65][fol - 65][2] + "\t");
						
						//maxmin phi
						out.append(plusplus[amino - 65][pre - 65][fol - 65][3] + "\t");
						out.append(plusplus[amino - 65][pre - 65][fol - 65][5] + "\t");
						
						//maxmin psi
						out.append(plusplus[amino - 65][pre - 65][fol - 65][4] + "\t");
						out.append(plusplus[amino - 65][pre - 65][fol - 65][6] + "\n");
						}else{
							out.append(999.0 + "\t"); //profile not exists
							out.append(999.0 + "\t" + 0.0 + "\t");
							out.append(360.0 + "\t" + -360.0 + "\t" );
							out.append(360.0 + "\t" + -360.0 + "\n" );
						}
					}
				}
			}
		out.close();
		} catch (IOException e) {
			System.out.println("cannot write angles ");
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
		DihedralAngles.calc("../GoBi_old/DSSP");
	}

}
