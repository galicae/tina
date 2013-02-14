package bioinfo.energy.potential.voronoi;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.OutputStreamWriter;

public class VoroPPWrap {
	
	private final String TMPDIR;
	private final String VOROBIN;

	public VoroPPWrap(String tmpdir, String vorobin){
		this.TMPDIR = tmpdir;
		this.VOROBIN = vorobin;
	}
	
	public VoroPPWrap(String vorobin){
		this.TMPDIR = "/tmp/";
		this.VOROBIN = vorobin;
	}
	
	/**
	 * CAVE: an extrenal tool will be used, it should be positioned where VOROBIN leads to
	 * CAVE: the extrenal tool must use files as input and output, so you should be sure that TMPDIR exists and is writeable
	 * decomposites the points of data and writes the results back to data
	 * @param data input points in value object of type VoronoiData
	 * @return data instance of VoronoiData, which is just the same as input data with some additional values 
	 */
	public void decomposite(VoronoiData data){
				
		String tmpfile = TMPDIR+data.getID()+".lst";
		File tmpfileF = new File(tmpfile);
		try{
			BufferedWriter bw = new BufferedWriter(new OutputStreamWriter(new FileOutputStream(tmpfile)));
			bw.append(data.toVoroPPString());
			bw.flush();
			bw.close();
		}catch(Exception e){
			e.printStackTrace();
		}
		
		Process pHull = null;;
		try{
			double[][] hull = data.calcHull();
			String sHull = hull[0][0]+" "+hull[0][1]+" "+hull[1][0]+" "+hull[1][1]+" "+hull[2][0]+" "+hull[2][1];
			//System.out.println("\twill process -> "+VOROBIN+" -c %i#%n#%f "+sHull+" "+tmpfile);
			pHull = Runtime.getRuntime().exec(VOROBIN+" -c %i#%n#%f "+sHull+" "+tmpfile);
			pHull.waitFor();
		
			String line;
	    	BufferedReader bri = new BufferedReader(new InputStreamReader(pHull.getInputStream()));
	    	BufferedReader bre = new BufferedReader(new InputStreamReader(pHull.getErrorStream()));
	    	//System.out.println("<error>");
	    	try{
		    	while ((line = bri.readLine()) != null) {
		    		System.err.println(line);
		    	}
		    	bri.close();
		    	while ((line = bre.readLine()) != null) {
		    		System.err.println(line);
		    	}
		    	bre.close();
		    	
	    	}catch(Exception ex){
	    		ex.printStackTrace();
	    	}
		} catch (Exception e){	
			e.printStackTrace();
		}
		//System.out.println("<error>");
		
		String resultfile = tmpfile+".vol";
		File resultfileF = new File(resultfile);
		if(resultfileF.exists()){
			try{
				BufferedReader br = new BufferedReader(new InputStreamReader(new FileInputStream(resultfile)));
				String line = null;
				int point;
				String[] comp = null;
				String[] neighborsS = null;
				String[] areasS = null;
				int[] neighbors = null;
				double[] areas = null;
				while((line = br.readLine()) != null){
					line = line.trim();
					comp = line.split("#");
					point = Integer.parseInt(comp[0].trim());
					neighborsS = comp[1].trim().split("\\s");
					areasS = comp[2].trim().split("\\s");
					neighbors = new int[neighborsS.length];
					areas = new double[areasS.length];
					for(int i = 0; i != neighbors.length; i++){
						neighbors[i] = Integer.parseInt(neighborsS[i]);
					}
					for(int i = 0; i != areas.length; i++){
						areas[i] = Double.parseDouble(areasS[i]);
					}
					data.addValues(point, neighbors, areas);
				}
				br.close();
			}catch(Exception e){
				e.printStackTrace();
			}
		}

	    tmpfileF.delete();
	    resultfileF.delete();
	}
	

}
