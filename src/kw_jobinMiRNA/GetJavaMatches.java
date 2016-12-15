package kw_jobinMiRNA;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashSet;

public class GetJavaMatches implements Runnable {
	public volatile boolean finished = false;
	public volatile Exception exception = null;

	private final File file;
	private final String outdir;
	private final ArrayList<String> keys;
	
	public final String id;
	public int numReads;
	public int numMatched;

	public GetJavaMatches(File f, String out, ArrayList<String> map) {
		file = f;
		outdir = out;
		keys = map;
		id = file.getName().replace(".adapterfiltered.fasta", "");
		numReads = 0;
		numMatched = 0;
	}

	@Override
	public void run() {
		try {
			HashSet<String> matched = new HashSet<String>();
			BufferedReader fa = new BufferedReader(new FileReader(file));
			BufferedWriter matchOut = new BufferedWriter(new FileWriter(new File(
					outdir + id + ".txt")));

			matchOut.write("referenceHeader\treadHeader\n");

			for(String line1 = fa.readLine(); line1 != null; line1 = fa.readLine()) {
				String seq = fa.readLine();
				String head = line1.replace(">", "");
				for(String r : keys) {
					if((seq.length() >= r.length() && seq.contains(r)) ||
							(seq.length() < r.length() && r.contains(seq))) {
						matched.add(head);
						matchOut.write(r + "\t" + head + "\n");
					}
				}
				numReads++;
			}
			fa.close();
			matchOut.close();	
			numMatched = matched.size();
		} catch (IOException e) {
			e.printStackTrace();
			this.exception = e;
		} finally {
			this.finished = true;
		}
	}
}
