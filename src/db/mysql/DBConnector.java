package db.mysql;

import java.sql.PreparedStatement;
import java.sql.ResultSet;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.List;

import bioinfo.proteins.AminoAcid;
import bioinfo.proteins.AminoAcidName;
import bioinfo.proteins.Atom;
import bioinfo.proteins.DSSPEntry;
import bioinfo.proteins.PDBEntry;
import bioinfo.proteins.SecStructEight;

public class DBConnector extends MysqlWrapper{

	private static final String pdb_table = "pdb";
	private static final String amino_table = "aminoacid";
	private static final String atom_table = "atom";
	private static final String dssp_table = "dssp";
	private static final String seq_table = "sequence";
	
	private static final String[] pdbfields = {"id","pdb_id","chainID", "chainIDNum","length"};
	private static final String[] dsspfields = {"id","ss","sa", "phi","psi","x_ca","y_ca","z_ca","aminoacid_id","pdb_id"};
	private static final String[] aafields = {"id","name","res_index","numberofAtom","pdb_id"};
	private static final String[] atomfields = {"id","type","x","y","z","aminoacid_id"};
	
	//pdb queries
	private static final String setPDBEntry = "insert into "+pdb_table+" ("+pdbfields[1]+","+pdbfields[2]+","+pdbfields[3]+","+pdbfields[4]+") values (?,?,?,?)";	
	private static final String getPDBById = "select * from "+pdb_table+" join aminoacid on pdb.id = aminoacid.pdb_id join atom on aminoacid.id = atom.aminoacid_id" +
			" where pdb.pdb_id = ? and pdb.chainID = ? and pdb.chainIDNum = ?";
	private static final String setAminoEntry = "insert into "+amino_table+" ("+aafields[1]+","+aafields[2]+","+aafields[3]+","+aafields[4]+") values (?,?,?,?)";
	private static final String setAtomEntry = "insert into "+atom_table+" ("+atomfields[1]+","+atomfields[2]+","+atomfields[3]+","+atomfields[4]+","+atomfields[5]+") values (?,?,?,?,?)";
	private static final String queryPDBExist = "select id from "+pdb_table+" where "+pdbfields[1]+" = ? and "+pdbfields[2]+" = ? and "+pdbfields[3]+" = ?";

	
	//dssp queries
	private static final String setDSSPEntry = "insert into "+dssp_table+" ("
			+ dsspfields[1] + "," + dsspfields[2] + "," + dsspfields[3] + ","
			+ dsspfields[4] + "," + dsspfields[5] + "," + dsspfields[6] + ","
			+ dsspfields[7] + "," + dsspfields[8] + "," + dsspfields[9] + ") values (?,?,?,?,?,?,?,?,?)";	
	private static final String getDSSPById = "select dssp.*,aminoacid.name,aminoacid.res_index from "+dssp_table+" join pdb on pdb.id = dssp.pdb_id join aminoacid on aminoacid.id = dssp.aminoacid_id" +
			" where pdb.pdb_id = ? and pdb.chainID = ? and pdb.chainIDNum = ?";
	private static final String queryDSSPExist = "select dssp.id from "+dssp_table+" join pdb on pdb.id = dssp.pdb_id where pdb."+pdbfields[1]+" = ? and "+pdbfields[2]+" = ? and "+pdbfields[3]+" = ?";
	
	
	//statements
	private PreparedStatement stmtGetPDBById = connection.createStatement(getPDBById);
	private PreparedStatement stmtGetDSSPById = connection.createStatement(getDSSPById);
	private PreparedStatement stmtPDBExist = connection.createStatement(queryPDBExist);
	private PreparedStatement stmtDSSPExist = connection.createStatement(queryDSSPExist);
	private PreparedStatement stmtLastIDAA = connection.createStatement("Select LAST_INSERT_ID() from "+amino_table);
	private PreparedStatement stmtLastIDPDB = connection.createStatement("Select LAST_INSERT_ID() from "+pdb_table);
	private PreparedStatement stmtSetPDBEntry = connection.createStatement(setPDBEntry);
	private PreparedStatement stmtSetAminoEntry = connection.createStatement(setAminoEntry);
	private PreparedStatement stmtSetAtomEntry = connection.createStatement(setAtomEntry);
	private PreparedStatement stmtSetDSSPEntry = connection.createStatement(setDSSPEntry);
	private PreparedStatement stmtGetAAIdByResIndex = connection.createStatement("select id from "+amino_table+" where res_index = ? and pdb_id = ?");
	
	public DBConnector(MysqlDBConnection connection) {
		super(connection);
	}
	
	@Override
	String getTablename() {
		return pdb_table;
	}

	@Override
	String[] getFields() {
		return pdbfields;
	}
	
	public PDBEntry getPDB(String id){
		List<Atom> atoms = new ArrayList<Atom>();
		List<AminoAcid> aminos = new ArrayList<AminoAcid>();
		
		String chainID = id.substring(4,5);
		int chainIDNum=Integer.parseInt(id.substring(5));
		String pdbid = id.substring(0, 4);
		
		try{
			stmtGetPDBById.setString(1, pdbid);
			stmtGetPDBById.setString(2, chainID);
			stmtGetPDBById.setInt(3, chainIDNum);
			ResultSet res = stmtGetPDBById.executeQuery();
			
			int atomquantity;
			double[] pos_temp;
			
			if(pdbExist(id) != -1){
				while(res.next()){
					//read out aminos and corresponding atoms
					atomquantity = res.getInt(aafields[3]);
					pos_temp = new double[3];
					pos_temp[0] = res.getDouble(atomfields[2]);
					pos_temp[1] = res.getDouble(atomfields[3]);
					pos_temp[2] = res.getDouble(atomfields[4]);
					atoms.add(new Atom(res.getString(atomfields[1]),pos_temp));
					
					for (int j = 0; j < atomquantity-1; j++) {
						res.next();
						pos_temp = new double[3];
						pos_temp[0] = res.getDouble(atomfields[2]);
						pos_temp[1] = res.getDouble(atomfields[3]);
						pos_temp[2] = res.getDouble(atomfields[4]);
						atoms.add(new Atom(res.getString(atomfields[1]),pos_temp));
					}
					aminos.add(new AminoAcid(AminoAcidName.getAAFromOLC(res.getString(aafields[1])),res.getInt(aafields[2]),atoms.toArray(new Atom[atomquantity])));
					atoms.clear();
				}
				return new PDBEntry(id,aminos.toArray(new AminoAcid[aminos.size()]));
			} else {
				return null;
			}
		}catch(SQLException e){
			e.printStackTrace();
			return null;
		}
	}
	
	public DSSPEntry getDSSP(String id){
		List<AminoAcidName> aminos = new ArrayList<AminoAcidName>();
		List<Integer> resIndex = new ArrayList<Integer>();
		List<SecStructEight> secStruct = new ArrayList<SecStructEight>();
		List<double[]> caTrace = new ArrayList<double[]>();
		List<Integer> sa = new ArrayList<Integer>();
		List<Double> phi = new ArrayList<Double>();
		List<Double> psi = new ArrayList<Double>();
		
		String pdbid = id.substring(0, 4);
		String chainID = id.substring(4,5);
		int chainIDNum=Integer.parseInt(id.substring(5));
		
		try{
			stmtGetDSSPById.setString(1, pdbid);
			stmtGetDSSPById.setString(2, chainID);
			stmtGetDSSPById.setInt(3, chainIDNum);
			ResultSet res = stmtGetDSSPById.executeQuery();

			if(pdbExist(id) != -1){		
				while(res.next()){
					//read out dssp
					aminos.add(AminoAcidName.getAAFromOLC(res.getString("aminoacid."+aafields[1])));
					resIndex.add(res.getInt("aminoacid."+aafields[2]));
					secStruct.add(SecStructEight.getSSFromChar(res.getString(dsspfields[1]).charAt(0)));
					caTrace.add(new double[] {res.getDouble(dsspfields[5]),res.getDouble(dsspfields[6]),res.getDouble(dsspfields[7])});
					sa.add(res.getInt(dsspfields[2]));
					phi.add(res.getDouble(dsspfields[3]));
					psi.add(res.getDouble(dsspfields[4]));
					
				}
				return new DSSPEntry(id,aminos,resIndex,secStruct,sa,phi,psi,caTrace);
			} else {
				return null;
			}
		}catch(SQLException e){
			e.printStackTrace();
			return null;
		}
	}
	
	public int pdbExist(String id){
		String pdbid = id.substring(0,4);
		String chain = id.substring(4,5);
		int chainIdNum = Integer.parseInt(id.substring(5));
		try{
			stmtPDBExist.setString(1, pdbid);
			stmtPDBExist.setString(2, chain);
			stmtPDBExist.setInt(3, chainIdNum);
			
			ResultSet res = stmtPDBExist.executeQuery();
			
			if(res.first()){
				return res.getInt(pdbfields[0]);
			}
			return -1;
		}catch(SQLException e){
			e.printStackTrace();
			return -1;
		}
	}
	
	public int dsspExist(String id){
		String pdbid = id.substring(0,4);
		String chain = id.substring(4,5);
		int chainIdNum = Integer.parseInt(id.substring(5));
		try{
			stmtDSSPExist.setString(1, pdbid);
			stmtDSSPExist.setString(2, chain);
			stmtDSSPExist.setInt(3, chainIdNum);
			
			ResultSet res = stmtDSSPExist.executeQuery();
			
			if(res.first()){
				return res.getInt(pdbfields[0]);
			}
			return -1;
		}catch(SQLException e){
			e.printStackTrace();
			return -1;
		}
	}
	
//	public String[] getIDs(){
//		Statement stmt = connection.createStatement();
//		try {
//			ResultSet res = stmt.executeQuery("select id from "+getTablename());
//			List<String> ids = new ArrayList<String>();
//			while(res.next()){
//				ids.add(res.getString(pdbfields[0]));
//			}
//			return ids.toArray(new String[ids.size()]);
//		} catch (SQLException e) {
//			e.printStackTrace();
//			return null;
//		}
//	}
		
	private int getLastIdAA(){
		int lastid;
		ResultSet res;
		try {
			res = stmtLastIDAA.executeQuery();
			res.first();
			lastid = res.getInt(1);
		} catch (SQLException e) {
			e.printStackTrace();
			return 0;
		}	
		return lastid;
	}
	
	private int getLastIdPDB(){
		int lastid;
		ResultSet res;
		try {
			res = stmtLastIDPDB.executeQuery();
			res.first();
			lastid = res.getInt(1);
		} catch (SQLException e) {
			e.printStackTrace();
			return 0;
		}	
		return lastid;
	}
	
	public boolean addPDBEntry(PDBEntry entry){		
		
		if(pdbExist(entry.getID()+entry.getChainID()+entry.getChainIDNum()) == -1){
			try{
				int lastpdbid;
				int lastaminoid;
				AminoAcid amino;
				Atom atom;
				double[] pos;
				
				//insert pdbentry
				stmtSetPDBEntry.setString(1,entry.getID());
				stmtSetPDBEntry.setString(2,String.valueOf(entry.getChainID()));
				stmtSetPDBEntry.setInt(3,entry.getChainIDNum());
				stmtSetPDBEntry.setInt(4, entry.length());
				stmtSetPDBEntry.execute();
				lastpdbid = getLastIdPDB();
				
				//insert aminos
				for (int i = 0; i < entry.length(); i++) {
					amino = entry.getAminoAcid(i);					
					stmtSetAminoEntry.setString(1, amino.toString());
					stmtSetAminoEntry.setInt(2, amino.getResIndex());
					stmtSetAminoEntry.setInt(3, amino.getAtomNumber());
					stmtSetAminoEntry.setInt(4, lastpdbid);
					stmtSetAminoEntry.execute();
					lastaminoid = getLastIdAA();
					
					//insert atoms
					for (int j = 0; j < amino.getAtomNumber(); j++) {
						atom = amino.getAtom(j);
						pos = atom.getPosition();
						stmtSetAtomEntry.setString(1, atom.getType().toString());
						stmtSetAtomEntry.setDouble(2, pos[0]);
						stmtSetAtomEntry.setDouble(3, pos[1]);
						stmtSetAtomEntry.setDouble(4, pos[2]);
						stmtSetAtomEntry.setInt(5, lastaminoid);
						stmtSetAtomEntry.addBatch();
					}
				}
				stmtSetAtomEntry.executeBatch();
				return true;		
			}catch(SQLException e){
				e.printStackTrace();
				return false;
			}
		} else {
			System.out.println("PDB already exists in DB");
			return false;
		}
	}
		
		public boolean addDSSPEntry(DSSPEntry entry){			
			int pdbDBID = pdbExist(entry.getID()+entry.getChainID()+entry.getChainIDNum());
			int aminoid;
			ResultSet res;
			
			if(pdbDBID != -1 && dsspExist(entry.getID()+entry.getChainID()+entry.getChainIDNum()) == -1){
				try{					
					for (int i = 0; i < entry.getLength(); i++) {
						stmtGetAAIdByResIndex.setInt(1, entry.getResIndex()[i]);
						stmtGetAAIdByResIndex.setInt(2, pdbDBID);
						res = stmtGetAAIdByResIndex.executeQuery();
						aminoid = -1;
						if(res.first()){
							aminoid = res.getInt(aafields[0]);
						}else{
							System.out.println(entry.getID()+entry.getChainID()+entry.getChainIDNum());
							System.out.println(entry.getResIndex()[i]);
							System.out.println("--------------------------");
						}
						
						stmtSetDSSPEntry.setString(1, entry.getSecondaryStructure()[i].toString());
						stmtSetDSSPEntry.setInt(2, entry.getAccesability()[i]);
						stmtSetDSSPEntry.setDouble(3, entry.getPhi()[i]);
						stmtSetDSSPEntry.setDouble(4, entry.getPsi()[i]);
						stmtSetDSSPEntry.setDouble(5, entry.getCaTrace()[i][0]);
						stmtSetDSSPEntry.setDouble(6, entry.getCaTrace()[i][1]);
						stmtSetDSSPEntry.setDouble(7, entry.getCaTrace()[i][2]);
						stmtSetDSSPEntry.setInt(8, aminoid);
						stmtSetDSSPEntry.setInt(9, pdbDBID);
						stmtSetDSSPEntry.addBatch();
					}
					stmtSetDSSPEntry.executeBatch();
					return true;
					
				}catch(SQLException e){
					e.printStackTrace();
					return false;
				}
			} else {
				System.out.println("Corresponding PDB entry doesn't exist or DSSP already imported.");
				return false;
			}
		}
}
