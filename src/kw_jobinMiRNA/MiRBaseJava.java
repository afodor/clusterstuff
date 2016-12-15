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

public class MiRBaseJava {
	public static String DIR = "/nobackup/afodor_research/kwinglee/jobin/microRNA/";
	public static String FQDIR = DIR + "adapterFiltered/";
	public static String OUTDIR = DIR + "miRBaseJava/";
	public static String MATREF = "/nobackup/afodor_research/kwinglee/mirbase_v21/mature.fa";
	public static String PINREF = "/nobackup/afodor_research/kwinglee/mirbase_v21/hairpin.fa";

	public static void main(String[] args) throws IOException {
		//set up reference based on input
		if(args.length != 1) {
			System.out.println("Input: 0 for mature or 1 for hairpin");
			System.exit(1);
		}
		String[] refs = new String[]{MATREF, PINREF};
		String[] refNames = new String[]{"mature", "hairpin"};
		String ref = refs[Integer.parseInt(args[0])];
		String refName = refNames[Integer.parseInt(args[0])];
		HashMap<String, String> refSeqs = new HashMap<String, String>();//map of sequence to key
		BufferedReader brRef = new BufferedReader(new FileReader(ref));
		String header = brRef.readLine();
		String line = brRef.readLine();
		String sequ = line;
		while(line != null) {
			if(line.startsWith(">")) {
				refSeqs.put(sequ, header.replace(">", ""));
				header = line;
				sequ = "";
				System.out.println(header);
				System.out.println(sequ);
			} else {
				sequ += line;
			}
			line = brRef.readLine();
		}
		refSeqs.put(sequ, header.replace(">", ""));
		brRef.close();
		ArrayList<String> keys = new ArrayList<String>(refSeqs.keySet());
		keys.remove(null);
		
		//for each read in each fasta file, see if contains an exact match for any database string
		BufferedWriter out = new BufferedWriter(new FileWriter(new File(
				OUTDIR + "MiRBaseJavaCounts_" + refName+ ".txt")));
		out.write("sampleID\tnumReads\tnumMatch\n");
		File[] files = new File(FQDIR).listFiles();
		for(File f : files) {
			if(f.getName().endsWith(".fasta")) {
				String id = f.getName().replace(".adapterfiltered.fasta", "");
				int numReads = 0;
				HashSet<String> matched = new HashSet<String>();
				BufferedReader fa = new BufferedReader(new FileReader(f));
				BufferedWriter matchOut = new BufferedWriter(new FileWriter(new File(
						OUTDIR + refName + "_v_" + id + ".txt")));
				matchOut.write("referenceHeader\treadHeader\n");
				for(String line1 = fa.readLine(); line1 != null; line1 = fa.readLine()) {
					String seq = fa.readLine();
					String head = line1.replace(">", "");
					for(String r : keys) {
						if((seq.length() >= r.length() && seq.contains(r)) ||
								(seq.length() < r.length() && r.contains(seq))) {
							matched.add(head);
							matchOut.write(refSeqs.get(r) + "\t" + head + "\n");
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
