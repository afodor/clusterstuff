/*
 * parse AbundantOTU+ results into OTU table
 * no chimeric sequences so don't need to filter
 * 5/17/16
 */
package kw_jobinDolphin;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;

public class ParseAbundantOTUResults {
	public static String DIR = "/nobackup/afodor_research/kwinglee/jobin/dolphin/abunOTU/";
	public static String QIIMEMAP = "/nobackup/afodor_research/kwinglee/jobin/dolphin/qiimeMap.txt";
	public static ArrayList<String> SAMPLES;
	
	public static void main(String[] args) throws IOException {
		//get list of samples
		SAMPLES = new ArrayList<String>();
		/*String[] files = new File("/nobackup/afodor_research/kwinglee/jobin/dolphin/stitched_reads").list();
		for(String f: files) {
			if(f.endsWith(".join.fasta")) {
				SAMPLES.add(f.replace(".join.fasta", "").replace("-", ""));//water samples are Water-Sample-# in stitched and WaterSample# in cluster file
			}
		}*/
		BufferedReader qiimeMap = new BufferedReader(new FileReader(new File(
				QIIMEMAP)));
		String line = qiimeMap.readLine();
		line = qiimeMap.readLine();
		while(line != null) {
			SAMPLES.add(line.split("\t")[0]);
			line = qiimeMap.readLine();
		}
		qiimeMap.close();
		Collections.sort(SAMPLES);
		
		//write header/set up output file
		BufferedWriter out = new BufferedWriter(new FileWriter(new File(
				DIR + "dolphinAbundantOTUtable.txt")));
		out.write("consensus");
		for(int i = 0; i < SAMPLES.size(); i++) {
			out.write("\t" + SAMPLES.get(i));
		}
		out.write("\n");
		
		/*
		 * clust file is set up as line starting with # indicating consensus
		 * number and number of sequences, followed by lines of which sequences clustered
		 * for each consensus, make map of sample to number of times sample in cluster
		 * then write results
		 */
		BufferedReader clust = new BufferedReader(new FileReader(new File(
				DIR + "dolphinAbundantOTU.clust")));
		HashMap<String, Integer> map = new HashMap<String, Integer>();
		String name = "";
		int numSeqs = 0;
		for(line = clust.readLine(); line != null; line=clust.readLine()) {
			if(line.startsWith("#")) {
				if(numSeqs > 0) {
					printCounts(out, name, numSeqs, map);
				}
				String[] sp = line.split(" ");
				name = "Consensus" + sp[1];
				numSeqs = Integer.parseInt(sp[3].replace("seq=", ""));
				map = new HashMap<String, Integer>();
			} else {
				String sample = line.trim().split(" ")[2].split("_")[0];
				if(map.containsKey(sample)) {
					map.put(sample, map.get(sample)+1);
				} else {
					map.put(sample, 1);
				}
			}
		}
		printCounts(out, name, numSeqs, map);
		
		clust.close();
		out.close();
	}
	
	public static void printCounts(BufferedWriter out, String cons, int numSeqs, 
			HashMap<String, Integer> map) throws IOException {
		out.write(cons);
		int check = 0;
		for(int i = 0; i < SAMPLES.size(); i++) {
			String samp = SAMPLES.get(i);
			if(map.containsKey(samp)) {
				out.write("\t" + map.get(samp));
				check += map.get(samp);
			} else {
				out.write("\t0");
			}
		}
		out.write("\n");
		out.flush();
		if(check != numSeqs) {
			System.err.println("Incorrect number of seqs: " + cons);
			System.err.println("  # seen: " + check);
			System.err.println("  # expected: " + numSeqs);
		}
	}

}
