/*
 * Parses alignment to cards protein homolog database
 * 6/2/16
 */
package kw_china_wgs;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;

public class ParseCardsAlignmentResults {
	public static String DIR = "/nobackup/afodor_research/kwinglee/china/wgs/";
	public static int NUM_SAMPS = 40;
	
	public static void main(String[] args) throws Exception {
		//get list of samples
		String[] fastas = new File(DIR + "fastas").list();
		String[] samples = new String[NUM_SAMPS];
		int check = 0;
		for(String f : fastas) {
			if(f.endsWith(".fa")) {
				samples[check] = f.replace("_1.fa", "");
				check++;
			}
		}
		if(check != NUM_SAMPS) {
			throw new Exception("Incorrect number of samples " + check);
		}
		
		//get counts (raw and non-human) for each gene
		HashMap<String, Integer[]> geneCounts = new HashMap<String, Integer[]>();
		HashMap<String, Integer[]> nonHumGeneCounts = new HashMap<String, Integer[]>();
		int[] numHits = new int[NUM_SAMPS];
		int[] numNonHumanHits = new int[NUM_SAMPS];
		
		//get number of hits (only focus on protein homolog results)
		String alignDir = DIR + "alignToCards/";
		for(int i = 0; i < samples.length; i++) {
			String samp = samples[i];
			//get set of reads that matched human genome
			HashSet<String> humanReads = getHumanReads(samp);
			//get reads that mapped to cards protein homolog
			BufferedReader sam = new BufferedReader(new FileReader(new File(
					alignDir + "pro_homolog_v_" + samp + ".mapped.sam")));
			//set up counts
			int numNotPrimary = 0;
			for(String read = sam.readLine(); read != null; read = sam.readLine()){
				if(!read.startsWith("@")) {
					String[] sp = read.split("\t");
					//only include read if primary alignment
					if(sp[1].equals("256")) {
						numNotPrimary++;
					} else {
						numHits[i]++;
						String gene = sp[2];

						if(!geneCounts.containsKey(gene)) {
							Integer[] counts = new Integer[NUM_SAMPS];
							Arrays.fill(counts, 0);
							geneCounts.put(gene, counts);
						}
						Integer[] counts = geneCounts.get(gene);
						counts[i]++;

						//reads that didn't map to human
						if(!humanReads.contains(sp[0])) {
							numNonHumanHits[i]++;

							if(!nonHumGeneCounts.containsKey(gene)) {
								Integer[] nonHumCounts = new Integer[NUM_SAMPS];
								Arrays.fill(nonHumCounts, 0);
								nonHumGeneCounts.put(gene, nonHumCounts);
							}
							Integer[] nonHumCounts = nonHumGeneCounts.get(gene);
							nonHumCounts[i]++;
						}
					}
				}
			}
			sam.close();
			System.out.println(samp + "\t" + numNotPrimary);
		}
		
		//write results
		BufferedWriter out = new BufferedWriter(new FileWriter(new File(
				alignDir + "pro_homolog_results.txt")));
		ArrayList<String> genes = new ArrayList<String>(geneCounts.keySet());
		Collections.sort(genes);
		
		//write header
		out.write("sample\tnumberReadsMapped\tnumberNonHumanReadsMapped");
		for(String g : genes) {
			out.write("\t" + g + "\t" + "nonhuman " + g);
		}
		out.write("\n");
		
		//write counts
		for(int i = 0; i < NUM_SAMPS; i++) {
			out.write(samples[i] + "\t" + numHits[i] + "\t" + numNonHumanHits[i]);
			for(String g : genes) {
				out.write("\t" + geneCounts.get(g)[i] + "\t");
				if(nonHumGeneCounts.containsKey(g)) {
					out.write(nonHumGeneCounts.get(g)[i]);
				} else {
					out.write(0);
				}
			}
			out.write("\n");
		}
		
		out.close();
		
	}

	/*
	 * for the given sample, return the set of reads that mapped to the human genome
	 */
	private static HashSet<String> getHumanReads(String sample) throws IOException {
		BufferedReader sam = new BufferedReader(new FileReader(new File(
				DIR + "alignToHG38/" + sample + "_1.hg38.mapped.sam")));
		HashSet<String> set = new HashSet<String>();
		String line = sam.readLine();
		while(line != null) {
			String[] sp = line.split("\t");
			set.add(sp[0]);
			line = sam.readLine();
		}
		sam.close();
		return(set);
	}
}
