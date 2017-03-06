/*
 * get counts for each small RNA in a given database with bowtie2 alignment
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
import java.util.HashSet;

public class IndivRNAbowtie {
	private static String DIR = "/nobackup/afodor_research/kwinglee/jobin/microRNA/";
	private static int NUM_SAMPS = 10;

	public static void main(String[] args) throws IOException {
		getCounts(DIR + "miRBaseBowtie/", "IndivMiRBaseHairpinBowtie");
		getCounts(DIR + "piRBaseBowtie/", "IndivPiRBaseBowtie");
	}

	private static void getCounts(String dir, String outFile) throws IOException {
		HashMap<String, ArrayList<Integer>> rnaCounts = new HashMap<String, ArrayList<Integer>>();//map of each rna to its count for each sample
		File[] files = new File(dir).listFiles();
		for(File f : files) {
			if((dir.contains("miR") && f.getName().endsWith("bowtie.hairpin.sam")) ||
					(dir.contains("piR") && f.getName().endsWith("piR.bowtie.sam"))) {
				HashMap<String, String> bestHit = new HashMap<String,String>();//map of read to the name of the best RNA match
				HashMap<String, Integer> bestQual = new HashMap<String, Integer>();//map of read to the best quality map
				BufferedReader res = new BufferedReader(new FileReader(f));
				HashSet<String> mapped = new HashSet<String>();//set of reads mapped (in case read mapped multiple times)
				for(String line = res.readLine(); line!=null; line=res.readLine()) {
					if(!line.startsWith("@")) {
						String[] sp = line.split("\\s+");
						if(!sp[1].equals("4")) { //read is mapped
							mapped.add(sp[0]);
							if(!sp[1].equals("256")) {//primary alignment
								int mapq = Integer.parseInt(sp[4]);
								String read = sp[0];
								if(!bestQual.containsKey(read)) {
									bestQual.put(read, mapq);
									bestHit.put(read, sp[2]);								
								} else {//update if better score
									//account for unavailable quality (255)
									int prev = bestQual.get(read);
									if(prev == 255 || 
											(prev < mapq && mapq != 255)){
										bestQual.put(read, mapq);
										bestHit.put(read, sp[2]);	
										if(f.getName().equals("Sample1.adapterfiltered.piR.bowtie.sam")) {
											System.out.println("read\t" + sp[2]);
										}
									}
								}
							}
						}
					}
				}
				res.close();
				
				//look what reads not picked up
				System.out.println(f.getName() + "\t" + mapped.size() + 
						"\t" + bestHit.size() + "\t" + bestQual.size());

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
		BufferedWriter out = new BufferedWriter(new FileWriter(new File(dir + outFile + ".txt")));
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
