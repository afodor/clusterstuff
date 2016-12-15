/*
 * Look for exact matches between reads and miRBase
 */
package kw_jobinMiRNA;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.util.ArrayList;
import java.util.HashMap;

public class PiRBaseJava {
	public static String DIR = "/nobackup/afodor_research/kwinglee/jobin/microRNA/";
	public static String FQDIR = DIR + "adapterFiltered/";
	public static String OUTDIR = DIR + "piRBaseJava/";
	public static String REF = "/nobackup/afodor_research/kwinglee/piRBase_v1.0/piR_mouse_v1.0.fa";
	
	public static void main(String[] args) throws Exception {
		File odir = new File(OUTDIR);
		if(!odir.exists()) {
			odir.mkdirs();
		}

		HashMap<String, String> refSeqs = new HashMap<String, String>();//map of sequence to key
		BufferedReader brRef = new BufferedReader(new FileReader(REF));
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

		//set up multithreaded
		ArrayList<GetJavaMatches> threadList = new ArrayList<GetJavaMatches>();


		//for each read in each fasta file, see if contains an exact match for any database string
		File[] files = new File(FQDIR).listFiles();
		for(File f : files) {
			if(f.getName().endsWith(".fasta")) {
				GetJavaMatches match = new GetJavaMatches(f, OUTDIR + "piRBase_v_", keys);
				new Thread(match).start();
				threadList.add(match);
			}
		}

		//wait to finish
		boolean wait = true;
		while(wait) {
			Thread.sleep(1000);
			wait = false;
			for(GetJavaMatches gjm : threadList) {
				if(!gjm.finished) {
					wait = true;
				}
			}
		}

		//write summary results
		BufferedWriter out = new BufferedWriter(new FileWriter(new File(
				OUTDIR + "PiRBaseJavaCounts.txt")));
		out.write("sampleID\tnumReads\tnumMatch\n");
		for(GetJavaMatches gjm : threadList) {
			out.write(gjm.id + "\t" + gjm.numReads + "\t" + gjm.numMatched + "\n");
		}
		out.close();

		//check nothing failed
		for(GetJavaMatches gjm : threadList) {
			if(gjm.exception != null) {
				throw new Exception(gjm.exception);
			}
		}	
	}
}

