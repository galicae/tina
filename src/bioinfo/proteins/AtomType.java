package bioinfo.proteins;

/**
 * 
 * @author gobi_4
 * 
 */
public enum AtomType {
	// U is for unknown
	U, C, N, O,
	CA, CB, CG, CD, CE, CZ, CH,
	NA, NB, NG, ND, NE, NZ, NH,
	OA, OB, OG, OD, OE, OZ, OH,
	
	CA1, CB1, CG1, CD1, CE1, CZ1, CH1,
	NA1, NB1, NG1, ND1, NE1, NZ1, NH1,
	OA1, OB1, OG1, OD1, OE1, OZ1, OH1,
	
	CA2, CB2, CG2, CD2, CE2, CZ2, CH2,
	NA2, NB2, NG2, ND2, NE2, NZ2, NH2,
	OA2, OB2, OG2, OD2, OE2, OZ2, OH2,
	
	CA3, CB3, CG3, CD3, CE3, CZ3, CH3,
	NA3, NB3, NG3, ND3, NE3, NZ3, NH3,
	OA3, OB3, OG3, OD3, OE3, OZ3, OH3,
	
	CA4, CB4, CG4, CD4, CE4, CZ4, CH4,
	NA4, NB4, NG4, ND4, NE4, NZ4, NH4,
	OA4, OB4, OG4, OD4, OE4, OZ4, OH4;	
	
	public static AtomType createFromString(String type){
		try{
			return AtomType.valueOf(type);
		} catch(Exception e){
			return AtomType.U;
		}
	}
}