package webservice.workers;

import java.io.BufferedInputStream;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;

/**
 * @author huberste
 * @lastchange 2013-02-26
 */
public class BaselineWorker extends Worker {

	private String ali;
	private String filter;

	private String result;

	public BaselineWorker(String jobFile) {
		super(jobFile);
	}

	@Override
	public void work() {
		// read WORKING job file
		readFile();

		try {
			String workingDir = "/home/h/huberste/gobi/abgabe";
			String arg1 = "-jar";
			String arg2 = "baseline.jar";
			String arg3 = "--ali";
			String arg4 = ali;
			String arg5 = "--filter";
			String arg6 = filter;
			String arg7 = "--ids";
			String arg8 = "/home/h/huberste/gobi/data/cathscop.ids";
			String arg9 = "-mss";
			String arg10 = "/home/h/huberste/gobi/data/cutoff05";
			String arg11 = "--pdb";
			String arg12 = "/home/h/huberste/gobi/data/pdb";
			String arg13 = "-o";
			String arg14 = DONE_FILE + ".pdb";

			ProcessBuilder pb = new ProcessBuilder("java", arg1, arg2, arg3,
					arg4, arg5, arg6, arg7, arg8, arg9, arg10, arg11, arg12,
					arg13, arg14);
			pb.directory(new File(workingDir));

			result = arg14;
			Process proc = pb.start();

			BufferedInputStream outstr = new BufferedInputStream(
					proc.getInputStream());
			byte[] buf = new byte[1024];
			int nr = outstr.read(buf);
			while (nr != -1) {
				// DO NOTHING!
			}
		} catch (Exception e) {
			e.printStackTrace();
		}

		// Write DONE job file
		writeResult();
	}

	@Override
	protected void readFile() {
		BufferedReader from = null;

		String line = null;
		try {
			from = new BufferedReader(new FileReader(JOB_FILE));
			while ((line = from.readLine()) != null) {
				if (line.startsWith("ALIGNMENT=")) {
					String alignment = line.substring(10);
					line = from.readLine();
					alignment += "\n"+line;
					line = from.readLine();
					alignment += "\n"+line + "\n";
					BufferedWriter bw = new BufferedWriter (new FileWriter(DONE_FILE + ".ali"));
					bw.write(alignment);
					bw.close();
					ali = DONE_FILE + ".ali";
				} else if (line.startsWith("FILTER=")) {
					filter = line.substring(7);
				}
			}
		} catch (IOException e) {
			System.err.println("Error while trying to read " + JOB_FILE + ".");
			e.printStackTrace();
		} finally {
			try {
				from.close();
			} catch (IOException e) {
				System.err
						.println("Error while trying close " + JOB_FILE + ".");
				e.printStackTrace();
			}
		}
	}

	@Override
	protected void writeResult() {
		BufferedReader from = null;
		BufferedWriter to = null;

		String line = null;
		try {
			from = new BufferedReader(new FileReader(JOB_FILE));
			to = new BufferedWriter(new FileWriter(DONE_FILE));
			while ((line = from.readLine()) != null) {
				to.write(line + "\n");
			}
			to.write("RESULT=\n");
			to.write(result);
		} catch (IOException e) {
			System.err.println("Error while trying to copy " + JOB_FILE
					+ " to " + DONE_FILE + ".");
			e.printStackTrace();
		} finally {
			try {
				if (from != null)
					from.close();
				if (to != null)
					to.close();
			} catch (IOException e) {
				System.err.println("Error while trying close FileStreams");
				e.printStackTrace();
			}
		}
		new File(JOB_FILE).delete();
	}

}
