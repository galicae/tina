package bioinfo.energy.potential.eval;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStreamWriter;
import java.io.Writer;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

import opt.CommandLine;
import opt.ParameterType;
import opt.Setting;
import opt.ValidationRule;

import bioinfo.energy.potential.AVoroPotential;
import bioinfo.energy.potential.DSSPSolventPotential;
import bioinfo.energy.potential.GridSolventPotential;
import bioinfo.energy.potential.IEnergy;
import bioinfo.energy.potential.SipplContactPotential;
import bioinfo.proteins.AminoAcidName;
import bioinfo.proteins.PDBEntry;
import bioinfo.proteins.PDBFileReader;
import bioinfo.proteins.structure.PDBSequenceShuffle;

public class Eval {
	
	/**
	 * generates evaluation test values to be plotted
	 * @param shuffle PDBSequenceShuffle instance with existing Distribution
	 * @param model to be shuffled numberOfShuffles times
	 * @param potential to be evaluated, should be prepared for preparation before calling it here
	 * this means, that prepare Sequence Scoring is called in here, all values that have to be set should be set here
	 * @param numberOfShuffles 
	 * @param output Writer stream to print the results to
	 */
	public static void nativeVsShuffles(PDBSequenceShuffle shuffle, PDBEntry model, IEnergy potential, int numberOfShuffles, Writer output){
		potential.prepareSequenceScoring(model);
		try{
			output.append("identifier\tenergy\n");
			output.append(model.getID()+"\t"+potential.getNativeScoring()+"\n");
		}catch(Exception e){
			e.printStackTrace();
		}
		AminoAcidName[] shuffled;
		for(int i = 0; i != numberOfShuffles; i++){
			shuffled = shuffle.shuffleSequence(model.length());
			try{
				output.append(model.getID()+"_shuffled"+i+"\t"+potential.getSequenceScoring(shuffled)+"\n");
			}catch(Exception e){
				e.printStackTrace();
			}
		}
		try{
			output.flush();
		}catch(Exception e){
			e.printStackTrace();
		}
	}
	
	public static void nativesVsShuffle(PDBSequenceShuffle shuffle, PDBFileReader reader, List<String> models, IEnergy potential, int numberOfShuffles, Writer output){
		try{
			output.append("identifier\tnativeenergy");
			for(int i = 0; i != numberOfShuffles; i++){
				output.append("\tshuffle"+i+"energy");
			}
			output.append("\n");
		}catch(Exception e){
			e.printStackTrace();
		}
		AminoAcidName[] shuffled;
		PDBEntry model;
		String mId;
		for(int m = 0; m != models.size(); m++){
			mId = models.get(m);
			model = reader.readFromFolderById(mId);
			potential.prepareSequenceScoring(model);
			try{
				output.append(model.getID()+"\t"+potential.getNativeScoring());
			}catch(Exception e){
				e.printStackTrace();
			}
			for(int i = 0; i != numberOfShuffles; i++){
				shuffled = shuffle.shuffleSequence(model.length());
				try{
					output.append("\t"+potential.getSequenceScoring(shuffled));
				}catch(Exception e){
					e.printStackTrace();
				}
			}
			try{
				output.append("\n");
			}catch(Exception e){
				e.printStackTrace();
			}
		}
		try{
			output.flush();
		}catch(Exception e){
			e.printStackTrace();
		}
	}
	
	public static void footprints(PDBFileReader reader, List<String> models, IEnergy potential, Writer output){
		double[][] data = new double[models.size()][];
		int maxlength = 0;
		int tmplength;
		String[] ids = new String[models.size()];
		PDBEntry model;
		
		for(int i = 0; i != models.size(); i++){
			ids[i] = models.get(i);
			model = reader.readFromFolderById(models.get(i));
			data[i] = potential.getAminoScores(model);
			if((tmplength = data[i].length) > maxlength){
				maxlength = tmplength;
			}
		}
		
		try{
			output.append(ids[0]+"");
			for(int i = 1; i != data.length; i++){
				output.append("\t"+ids[i]);
			}
			output.append("\n");
			
			for(int l = 0; l != maxlength; l++){
				if(data[0].length > l){
					output.append(data[0][l]+"");
				}
				for(int i = 1; i != data.length; i++){
					if(data[i].length > l){
						output.append("\t"+data[i][l]);
					}else{
						output.append("\t");
					}
				}
				output.append("\n");
			}
			
			output.flush();
		}catch(Exception e){
			e.printStackTrace();
		}
		
	}
	
	public static void footprintVsShuffles(PDBSequenceShuffle shuffle, PDBEntry model, IEnergy potential, int numberOfShuffles, Writer output){
		potential.prepareSequenceScoring(model);
		double[][] data = new double[numberOfShuffles+1][model.length()];
		String[] ids = new String[numberOfShuffles+1];
		data[0] = potential.getNativeAminoScoring();
		ids[0] = model.getID()+model.getChainID()+String.format("%02d", model.getChainIDNum());
		
		AminoAcidName[] shuffled;
		for(int i = 0; i != numberOfShuffles; i++){
			shuffled = shuffle.shuffleSequence(model.length());
			data[i+1] = potential.getSequenceAminoScoring(shuffled);
			ids[i+1] = ids[0]+"_shuffle"+i;
		}
		
		try{
			output.append(ids[0]);
			for(int j = 1; j != numberOfShuffles+1; j++){
				output.append("\t"+ids[j]);
			}
			output.append("\n");
			
			for(int i = 0; i != model.length(); i++){
				output.append(data[0][i]+"");
				for(int j = 1; j != numberOfShuffles+1; j++){
					output.append("\t"+data[j][i]);
				}
				output.append("\n");
			}
			output.flush();
		}catch(Exception e){
			e.printStackTrace();
		}
	}
	

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		
		ValidationRule potentialRule = new ValidationRule() {
			@Override
			public boolean validate(String value) {
				List<String> potentialtypes = new ArrayList<String>();
				potentialtypes.add("sippl");
				potentialtypes.add("dsspsol");
				potentialtypes.add("gridsol");
				if(potentialtypes.contains(value)){
					return true;
				}
				return false;
			}
		};
		List<ValidationRule> potentialRules = new ArrayList<ValidationRule>();
		potentialRules.add(potentialRule);
		
		ValidationRule benchmarkRule = new ValidationRule() {
			@Override
			public boolean validate(String value) {
				List<String> potentialtypes = new ArrayList<String>();
				potentialtypes.add("vsrandom");
				potentialtypes.add("footprint");
				if(potentialtypes.contains(value)){
					return true;
				}
				return false;
			}
		};
		List<ValidationRule> benchmarkRules = new ArrayList<ValidationRule>();
		benchmarkRules.add(benchmarkRule);
		
		ValidationRule shuffleRule = new ValidationRule() {
			@Override
			public boolean validate(String value) {
				List<String> potentialtypes = new ArrayList<String>();
				potentialtypes.add("random");
				potentialtypes.add("background");
				if(potentialtypes.contains(value)){
					return true;
				}
				return false;
			}
		};
		List<ValidationRule> shuffleRules = new ArrayList<ValidationRule>();
		shuffleRules.add(shuffleRule);
		
		List<Setting> settings = new ArrayList<Setting>();
		settings.add(new Setting("-data", "data", true, null, ParameterType.String, false, null));
		settings.add(new Setting("-potential", "potential",true,null,ParameterType.String, false,null));
		settings.add(new Setting("-vorobin", "vorobin",true,null,ParameterType.String,false,null));
		settings.add(new Setting("-minContact","minContact",false,null,ParameterType.Double,false,null));
		settings.add(new Setting("-gridExtend","gridExtend",false,null,ParameterType.Double,false,null));
		settings.add(new Setting("-gridDensity","gridDensity",false,null,ParameterType.Double,false,null));
		settings.add(new Setting("-gridClash","gridClash",false,null,ParameterType.Double,false,null));
		settings.add(new Setting("-shuffles","shuffles",true,2,ParameterType.Integer,false,null));
		settings.add(new Setting("-out","out",false,null,ParameterType.String,false,null));
		settings.add(new Setting("-model","model",false,null,ParameterType.String,false,null));
		settings.add(new Setting("-pottype","pottype",true,null,ParameterType.String,true,potentialRules));
		settings.add(new Setting("-benchtype","benchtype",true,null,ParameterType.String,true,benchmarkRules));
		settings.add(new Setting("-shuffletype","shuffletype",true,null,ParameterType.String,true,shuffleRules));
		settings.add(new Setting("-h","help",false,null,ParameterType.Boolean,false,null));
		
		CommandLine command = null;
		try{
			command = new CommandLine(settings, args);
		}catch(Exception e){
			e.printStackTrace();
			System.out.println(Eval.produceHelp());
			return;
		}
		HashMap<String, Object> coms = command.getParsedComs();
		
		if(coms.get("help") != null){
			System.out.println(Eval.produceHelp());
			return;
		}
		
		String dataFolder = (String)coms.get("data");
		String potentialFile = (String)coms.get("potential");
		String vorobin = (String)coms.get("vorobin");
		
		String shuffletype = (String)coms.get("shuffletype");
		boolean voroParamSet = false;
		double minContact = 0.0d;
		double gridExtend = 0.0d;
		double gridDensity = 0.0d;
		double gridClash = 0.0d;
		
		if(coms.keySet().contains("minContact") && coms.keySet().contains("gridExtend") &&
				coms.keySet().contains("gridDensity") && coms.keySet().contains("gridClash")){
			minContact = (Double)coms.get("minContact");
			gridExtend = (Double)coms.get("gridExtend");
			gridDensity = (Double)coms.get("gridDensity");
			gridClash = (Double)coms.get("gridClash");
			voroParamSet = true;
		}
	
		int numberOfShuffles = (Integer)coms.get("shuffles");
		
		Writer output = null;
		PDBSequenceShuffle shuffle;
		if(coms.get("out") != null){
			String outfile = (String)coms.get("out");
			try {
				output = new BufferedWriter(new OutputStreamWriter(new FileOutputStream(outfile)));
			} catch (FileNotFoundException e) {
				e.printStackTrace();
			}
		}else{
			output = new BufferedWriter(new OutputStreamWriter(System.out));
		}
				
		IEnergy potential = null;;
		String pottype = (String)coms.get("pottype");
		if(pottype.equals("sippl")){
			potential = new SipplContactPotential();
		}else if (pottype.equals("dsspsol")){
			if(voroParamSet){
				potential = new DSSPSolventPotential(vorobin, minContact, gridExtend, gridDensity, gridClash);
				((AVoroPotential)potential).setVoroParam(vorobin, minContact, gridExtend, gridDensity, gridClash);
			}else{
				System.out.println("Please set voro params!");
				System.out.println(Eval.produceHelp());
			}
		}else if (pottype.equals("gridsol")){
			if(voroParamSet){
				potential = new GridSolventPotential(vorobin, minContact, gridExtend, gridDensity, gridClash);
				((AVoroPotential)potential).setVoroParam(vorobin, minContact, gridExtend, gridDensity, gridClash);
			}else{
				System.out.println("Please set voro params!");
				System.out.println(Eval.produceHelp());
			}
		}else{
			System.err.println("pottype unknown!");
		}
 
		potential.readFromFile(potentialFile);
		PDBFileReader reader = new PDBFileReader(dataFolder);
		String modelId = null;
		
		if(((String)coms.get("benchtype")).equals("vsrandom")){
			if((modelId = (String)coms.get("model")) == null){
				List<String> models = new ArrayList<String>();
				File file = new File(dataFolder);
				for (File f : file.listFiles()) {
					models.add(f.getName().substring(0, 7));
				}
				if(shuffletype.equals("random")){
					shuffle = new PDBSequenceShuffle();
				}else if(shuffletype.equals("background")){
					shuffle = new PDBSequenceShuffle(reader,models);
				}else{
					shuffle = null;
				}
				Eval.nativesVsShuffle(shuffle, reader, models, potential, numberOfShuffles, output);
			}else{
				PDBEntry model = reader.readFromFolderById(modelId);
				if(shuffletype.equals("random")){
					shuffle = new PDBSequenceShuffle();
				}else if(shuffletype.equals("background")){
					shuffle = new PDBSequenceShuffle(model);
				}else{
					shuffle = null;
				}
				Eval.nativeVsShuffles(shuffle, model, potential, numberOfShuffles, output);
			}
		}else if(((String)coms.get("benchtype")).equals("footprint")){
			if((modelId = (String)coms.get("model")) == null){
				List<String> models = new ArrayList<String>();
				File file = new File(dataFolder);
				for (File f : file.listFiles()) {
					models.add(f.getName().substring(0, 7));
				}
				if(shuffletype.equals("random")){
					shuffle = new PDBSequenceShuffle();
				}else if(shuffletype.equals("background")){
					shuffle = new PDBSequenceShuffle(reader,models);
				}else{
					shuffle = null;
				}
				Eval.footprints(reader, models, potential, output);
			}else{
				PDBEntry model = reader.readFromFolderById(modelId);
				if(shuffletype.equals("random")){
					shuffle = new PDBSequenceShuffle();
				}else if(shuffletype.equals("background")){
					shuffle = new PDBSequenceShuffle(model);
				}else{
					shuffle = null;
				}
				Eval.footprintVsShuffles(shuffle, model, potential, numberOfShuffles, output);
			}
		}else{
			System.err.println("benchtype unknown!");

		}
		
		try {
			output.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
		
		
		
	}
	
	public static String produceHelp(){
		return "This tool can produce several values for evaluation of energy potentials. \n" +
				"Note the following parameters:\n" +
				"\tNAME              REQUIRED     DESCRIPTION\n" +
				"\tdata              YES          location of datas in pdb format\n" +
				"\tpotential         YES          potential .pot file\n" +
				"\tpottype           YES          one of sippl, dsspsol, gridsol\n"+
				"\tbenchtype         YES          one of vsrandom, footprint\n" +
				"\tmodel             NO           identifier of model to evaluate\n" +
				"\tvorobin           YES          location of voro++ binary\n" +
				"\tout               NO           output filename, else stdout\n" +
				"\tminContact        YES          value for minContact\n" +
				"\tgridExtend        YES          extended hull of voro decomp, CAVE: allways higher then clash value!\n" +
				"\tgridDensity       YES          density of solvent points\n" +
				"\tgridClash         YES          clash value for peptide-solvent points\n" +
				"\tshuffles          NO           number of shuffles for evaluation (std 2), sometimes not needed eg \n" +
				"                                      for footprints without any model defined\n" +
				"\tshuffletype       YES          one of random or background\n" +
				"\th                 NO           produces this help\n\n";
	}

}
