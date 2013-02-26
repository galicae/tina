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
public class MusterWorker extends Worker {

	private String template;
	private String target;

	private String result;

	public MusterWorker(String jobFile) {
		super(jobFile);
	}

	@Override
	public void work() {
		// read WORKING job file
		readFile();

		try {
			String workingDir = "/home/h/huberste/gobi/abgabe";
			String arg1 = "-jar";
			String arg2 = "muster.jar";
			String arg3 = template;
			String arg4 = target;

			ProcessBuilder pb = new ProcessBuilder("java", arg1, arg2, arg3,
					arg4);
			pb.directory(new File(workingDir));

			result = "";
			Process proc = pb.start();

			BufferedInputStream outstr = new BufferedInputStream(
					proc.getInputStream());
			byte[] buf = new byte[1024];
			int nr = outstr.read(buf);
			while (nr != -1) {
				for (int i = 0; i < nr; i++) {
					result += (char) buf[i];
				}
				nr = outstr.read(buf);
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
				if (line.startsWith("TEMPLATE_ID=")) {
					template = line.substring(12);
				} else if (line.startsWith("TARGET_ID=")) {
					target = line.substring(10);
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
