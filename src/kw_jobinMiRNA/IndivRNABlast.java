/*
 * get counts for each small RNA in a given database
 */
package kw_jobinMiRNA;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;

public class IndivRNABlast {
	private static String DIR = "/nobackup/afodor_research/kwinglee/jobin/microRNA/";
	private static int NUM_SAMPS = 10;

	public static double PID = 99;//minimum percent ID allowed
	public static int MISMATCH = 1;//maximum number of mismatches
	
	public static void main(String[] args) throws IOException {
		getCounts(DIR + "miRBaseBlast/", "IndivMiRBaseHairpinBlast");
		getCounts(DIR + "piRBaseBlast/", "IndivPiRBaseBlast");
	}

	public static void getCounts(String dir, String outFile) throws IOException {
		HashMap<String, ArrayList<Integer>> rnaCounts = new HashMap<String, ArrayList<Integer>>();//map of each rna to its count for each sample
		File[] files = new File(dir).listFiles();
		for(File f : files) {
			if((dir.contains("miR") && f.getName().endsWith("blast.hairpin.txt")) ||
					(dir.contains("piR") && f.getName().endsWith("piR.blast.txt"))) {
				//for each read, get the best matching small RNA
				HashMap<String, Double> bestBitScore = new HashMap<String, Double>();//map of read to best bitScore
				HashMap<String, String> bestHit = new HashMap<String,String>();//map of read to the name of the best RNA match
				BufferedReader res = new BufferedReader(new FileReader(f));
				for(String line = res.readLine(); line!=null; line=res.readLine()) {
					String[] sp = line.split("\t");
					if(sp.length > 5) {
						String read = sp[0];
						Double bit = Double.parseDouble(sp[11]);
						//qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore
						if(Double.parseDouble(sp[2]) > PID &&
								Integer.parseInt(sp[4]) <= MISMATCH) {
							if(!bestBitScore.containsKey(read) || 
									bestBitScore.get(read) < bit) { //update if read is not in map or if match has better bit score
								bestBitScore.put(read, bit);
								bestHit.put(read, sp[1]);								
							}
						}
					} 
				}
				res.close();
				
				//compile into counts for each sample
				int sample = Integer.parseInt(
						f.getName().split("\\.")[0].replace("Sample", ""))-1;//-1 for 0 based
				for(String read : bestHit.keySet()) {
					String rna = bestHit.get(read);
					if(!rnaCounts.containsKey(rna)) {
						ArrayList<Integer> counts = new ArrayList<Integer>(NUM_SAMPS);
						for(int c = 0; c < NUM_SAMPS; c++) {
							if(c == sample) {
								counts.add(1);
							} else {
								counts.add(0);
							}
						}
						rnaCounts.put(rna, counts);
					} else {
						ArrayList<Integer> counts = rnaCounts.get(rna);
						counts.set(sample, counts.get(sample)+1);
					}
				}
			}
		}
		
		//write results
		ArrayList<String> keys = new ArrayList<String>();
		keys.addAll(rnaCounts.keySet());
		Collections.sort(keys);
		BufferedWriter out = new BufferedWriter(new FileWriter(new File(dir + outFile)));
		//header
		out.write("Sample");
		for(String k : keys) {
			out.write("\t" + k);
		}
		out.write("\n");
		//counts
		for(int s = 0; s < NUM_SAMPS; s++) {
			out.write("Sample" + (s+1));
			for(String k : keys) {
				out.write("\t" + rnaCounts.get(k).get(s));
			}
			out.write("\n");
		}
		out.close();
	}
}
