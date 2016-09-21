/*
 * Parses alignment to virulence databases
 * 9/20/16
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

public class ParseVirulenceAlignmentResults {
	public static String DIR = "/nobackup/afodor_research/kwinglee/china/wgs/";
	public static int NUM_SAMPS = 40;

	public static void main(String[] args) throws Exception {
		//get list of samples
		String fastaDir = DIR + "fastas/";
		String[] fastas = new File(fastaDir).list();
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

		//get number of reads in each file
		double[] totReads = new double[NUM_SAMPS];
		for(int i = 0; i < NUM_SAMPS; i++) {
			String sample = samples[i];
			BufferedReader fa = new BufferedReader(new FileReader(new File(
					fastaDir + sample + "_1.fa")));
			for(String line = fa.readLine(); line != null; line = fa.readLine()) {
				if(line.startsWith(">")) {
					totReads[i]++;
				}
			}
			fa.close();
		}
		
		//get number of hits (only focus on protein homolog results)
		String alignDir = DIR + "alignToVirulenceDB/";
		String[] databases = new String[]{"VFDBfull", "VFDBcore", "MvirDB"};
		for(String db : databases) {
			System.out.println(db);
			
			//get counts (all and non-human) for each gene
			double[] totNonHumanReads = new double[NUM_SAMPS];
			HashMap<String, Integer[]> geneCounts = new HashMap<String, Integer[]>();
			HashMap<String, Integer[]> nonHumGeneCounts = new HashMap<String, Integer[]>();
			int[] numHits = new int[NUM_SAMPS];
			int[] numNonHumanHits = new int[NUM_SAMPS];

			for(int i = 0; i < samples.length; i++) {
				String samp = samples[i];
				//get set of reads that matched human genome
				HashSet<String> humanReads = getHumanReads(samp);
				totNonHumanReads[i] = totReads[i] - humanReads.size();
				//get reads that mapped to cards protein homolog
				BufferedReader sam = new BufferedReader(new FileReader(new File(
						alignDir + db + "_v_" + samp + ".mapped.sam")));
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
					alignDir + db + "_results.txt")));
			ArrayList<String> genes = new ArrayList<String>(geneCounts.keySet());
			Collections.sort(genes);

			//write header
			out.write("sample\tnumberReads\tnumberNonHumanReads\tproportionReadsMapped\tproportionReadsMappedNonHuman");
			for(String g : genes) {
				out.write("\t" + g + "\t" + "nonhuman " + g);
				//out.write("\t" + g);
			}
			out.write("\n");

			//write counts
			for(int i = 0; i < NUM_SAMPS; i++) {
				out.write(samples[i] + "\t" + totReads[i] + "\t" + totNonHumanReads[i] +
						"\t" + (numHits[i] / totReads[i]) + "\t" + 
						(numNonHumanHits[i] / totNonHumanReads[i]));
				if(numHits[i] != numNonHumanHits[i]) {
					System.err.println("human hits mapped in " + samples[i] + 
							" " + (numHits[i] - numNonHumanHits[i]));
				}
				for(String g : genes) {
					//out.write("\t" + geneCounts.get(g)[i]);
					out.write("\t" + (geneCounts.get(g)[i] / totReads[i]) + "\t");
					if(nonHumGeneCounts.containsKey(g)) {
						out.write(String.valueOf(nonHumGeneCounts.get(g)[i] / totNonHumanReads[i]));
					} else {
						out.write("0");
					}
				}
				out.write("\n");
			}

			out.close();
			System.out.println();
		}
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
