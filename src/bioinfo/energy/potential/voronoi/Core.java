package bioinfo.energy.potential.voronoi;

import java.util.HashMap;
import java.util.List;

import bioinfo.proteins.AminoAcidName;
import bioinfo.proteins.SecStructThree;

public class Core {
	
	private String peptideId;
	private SecStructThree type;
	private int sequentialNumber;
	private int length;
	private int firstResidue;
	private int lastResidue;
	private int contactCount;
	private HashMap<Integer,Integer> coreContacts; //contact-core-number , number of sequence contacts
	private AminoAcidName[] sequence;
	
	public Core(String peptideId, SecStructThree type, int sequentialNumber, int length, int firstResidue, int lastResidue, int contactCount, HashMap<Integer,Integer> coreContacts, AminoAcidName[] sequence) {
		super();
		this.peptideId = peptideId;
		this.type = type;
		this.sequentialNumber = sequentialNumber;
		this.length = length;
		this.firstResidue = firstResidue;
		this.lastResidue = lastResidue;
		this.contactCount = contactCount;
		this.coreContacts = coreContacts;
		this.sequence = sequence;
	}
	
	public String getPeptideId() {
		return peptideId;
	}
	public void setPeptideId(String peptideId) {
		this.peptideId = peptideId;
	}
	public SecStructThree getType() {
		return type;
	}
	public void setType(SecStructThree type) {
		this.type = type;
	}
	public int getSequentialNumber() {
		return sequentialNumber;
	}
	public void setSequentialNumber(int sequentialNumber) {
		this.sequentialNumber = sequentialNumber;
	}
	public int getLength() {
		return length;
	}
	public void setLength(int length) {
		this.length = length;
	}
	public int getFirstResidue() {
		return firstResidue;
	}
	public void setFirstResidue(int firstResidue) {
		this.firstResidue = firstResidue;
	}
	public int getLastResidue() {
		return lastResidue;
	}
	public void setLastResidue(int lastResidue) {
		this.lastResidue = lastResidue;
	}
	public int getContactCount() {
		return contactCount;
	}
	public void setContactCount(int contactCount) {
		this.contactCount = contactCount;
	}
	public HashMap<Integer,Integer> getCoreContacts() {
		return coreContacts;
	}
	public void setCoreContacts(HashMap<Integer,Integer> coreContacts) {
		this.coreContacts = coreContacts;
	}
	public AminoAcidName[] getSequence() {
		return sequence;
	}
	public void setSequence(AminoAcidName[] sequence) {
		this.sequence = sequence;
	}
	
	
	
}
