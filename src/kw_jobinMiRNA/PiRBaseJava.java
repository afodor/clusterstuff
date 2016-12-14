/*
 * Look for exact matches between reads and miRBase
 */
package kw_jobinMiRNA;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;

public class PiRBaseJava {
	public static String DIR = "/nobackup/afodor_research/kwinglee/jobin/microRNA/";
	public static String FQDIR = DIR + "adapterFiltered/";
	public static String OUTDIR = DIR + "piRBaseJava/";
	public static String REF = "/nobackup/afodor_research/kwinglee/piRBase_v1.0/piRbaseMouseBowtie";
	
	public static void main(String[] args) throws IOException {
		File odir = new File(OUTDIR);
		if(!odir.exists()) {
			odir.mkdirs();
		}
		
		HashMap<String, String> refSeqs = new HashMap<String, String>();//map of sequence to key
		BufferedReader brRef = new BufferedReader(new FileReader(REF));
		for(String line1 = brRef.readLine(); line1 != null; line1 = brRef.readLine()) {
			String line2 = brRef.readLine();
			refSeqs.put(line2, line1.replace(">", ""));
		}
		brRef.close();
		ArrayList<String> keys = new ArrayList<String>(refSeqs.keySet());
		keys.remove(null);
		
		//for each read in each fasta file, see if contains an exact match for any database string
		BufferedWriter out = new BufferedWriter(new FileWriter(new File(
				OUTDIR + "PiRBaseJavaCounts.txt")));
		out.write("sampleID\tnumReads\tnumMatch\n");
		File[] files = new File(FQDIR).listFiles();
		for(File f : files) {
			if(f.getName().endsWith(".fasta")) {
				String id = f.getName().replace(".adapterfiltered.fasta", "");
				int numReads = 0;
				HashSet<String> matched = new HashSet<String>();
				BufferedReader fa = new BufferedReader(new FileReader(f));
				BufferedWriter matchOut = new BufferedWriter(new FileWriter(new File(
						OUTDIR + "piRBase_v_" + id + ".txt")));
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
				out.write(id + '\t' + numReads + "\t" + matched.size() + "\n");
			}
		}
		out.close();		
	}
}
