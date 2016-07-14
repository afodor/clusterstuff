/*
 * parse the results of blasting chs genomes against cards protein homolog database
 * using the scaffold results
 */
package kw_cre;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Set;

public class ParseCardsBlastResultsScaffolds {
	public static int NUM_SAMPLES = 76;
	public static String DIR = "/nobackup/afodor_research/kwinglee/cre/chs_v_cards/";
	public static double PID_CUT = 80;//cutoff for percent identity
	public static double LEN_CUT = 90;//cutoff for length
	//public static int BASES_APART = 100;//minimum number of bases separating hits

	public static void main(String[] args) throws Exception {
		//get list of genomes
		String[] files = new File("/nobackup/afodor_research/af_broad/carolina/").list();
		String[] samples = new String[NUM_SAMPLES];
		int index = 0;
		for(String f : files) {
			if(f.endsWith(".genes.gtf")) {
				samples[index] = f.replace(".genes.gtf", "");
				index++;
			}
		}
		if(index != NUM_SAMPLES) {
			throw new Exception("Incorrect number of samples " + index);
		}

		//get list of cards genes
		//for each gene in cards database, get its length
		BufferedReader db = new BufferedReader(new FileReader(new File(
				"/users/kwinglee/card/nucleotide_fasta.protein_homolog.fasta")));
		HashMap<String, Integer> cardsLengths = new HashMap<String, Integer>();
		String gene = "";
		int length = 0;
		for(String line = db.readLine(); line != null; line = db.readLine()) {
			if(line.startsWith(">")) {
				if(!gene.equals("")) {
					cardsLengths.put(gene, length);
				}
				gene = line.split(" ")[0].replace(">", "");
				length = 0;
			} else {
				length += line.length();
			}
		}
		db.close();
		cardsLengths.put(gene, length);

		//for scaffold blasts, for each genome get set of regions that mapped to each cards gene
		HashMap<String, ArrayList<Set<String>>> hits = new HashMap<String, ArrayList<Set<String>>>();
		HashMap<String, ArrayList<Set<String>>> shortNameHits = new HashMap<String, ArrayList<Set<String>>>();
		for(int i = 0; i < samples.length; i++) {
			String samp = samples[i];
			int numHits = 0;
			BufferedReader br = new BufferedReader(new FileReader(new File(
					DIR + samp + ".scaffolds_v_cardsProHomolog")));
			for(String line = br.readLine(); line != null; line = br.readLine()) {
				if(!line.startsWith("#")) {
					String[] sp = line.split("\t");
					String cards = sp[1];
					double pid = Double.parseDouble(sp[2]);//percent identity
					double len = 100.0 * (Double.parseDouble(sp[9]) - 
							Double.parseDouble(sp[8])) / cardsLengths.get(cards);//cards length
					if(pid > PID_CUT && len > LEN_CUT) {
						String region = "[" + sp[0] + "," + sp[6] + "," + sp[7] +
								"," + sp[11].trim() + "]";//region is [scaffold, start, stop, bit score]
						numHits++;
						if(!hits.containsKey(cards)) {
							ArrayList<Set<String>> sets = new ArrayList<Set<String>>(NUM_SAMPLES);
							for(int j = 0; j < NUM_SAMPLES; j++) {
								sets.add(null);
							}
							hits.put(cards, sets);
						}
						ArrayList<Set<String>> sets = hits.get(cards);
						if(sets.get(i) == null) {
							sets.set(i,new HashSet<String>());
						}
						//check does not overlap any regions already found
						sets.set(i, checkOverlap(region, sets.get(i)));


						//also add to short name
						String[] sp2 = cards.split("\\|");
						String sName = sp2[sp2.length-1].split("-")[0];
						if(!shortNameHits.containsKey(sName)) {
							sets = new ArrayList<Set<String>>(NUM_SAMPLES);
							for(int j = 0; j < NUM_SAMPLES; j++) {
								sets.add(null);
							}
							shortNameHits.put(sName, sets);
						}
						sets = shortNameHits.get(sName);
						if(sets.get(i) == null) {
							sets.set(i,new HashSet<String>());
						}
						//check does not overlap any regions already found
						sets.set(i, checkOverlap(region, sets.get(i)));
					}
				}
			}
			System.out.println(samp + " " + numHits + " gene hits");
			//genomeHits.put(samp, new Integer(numHits));
			br.close();
		}

		//write results
		BufferedWriter out = new BufferedWriter(new FileWriter(new File(
				DIR + "blastCardsProHomologScaffoldTables_pid" + PID_CUT 
				+ "_len" + LEN_CUT + ".txt")));
		//header
		out.write("CARDSgene");
		for(int i = 0; i < samples.length; i++) {
			out.write("\t" + samples[i]);
		}
		out.write("\n");
		ArrayList<String> keys = new ArrayList<String>(hits.keySet());
		Collections.sort(keys);
		for(String k : keys) {
			out.write(k);
			ArrayList<Set<String>> matches = hits.get(k);
			for(int gen = 0; gen < matches.size(); gen++) {
				if(matches.get(gen) == null || matches.get(gen).size() == 0) {
					out.write("\tNA");
				} else {
					ArrayList<String> genes = new ArrayList<String>(matches.get(gen));
					out.write("\t" + genes.get(0));
					for(int g = 1; g < genes.size(); g++) {
						out.write(";" + genes.get(g));
					}
				}
			}
			out.write("\n");
		}
		out.close();

		//write collapse hits by gene type
		out = new BufferedWriter(new FileWriter(new File(
				DIR + "blastCardsProHomologScaffoldTables_pid" + PID_CUT 
				+ "_len" + LEN_CUT + "_collapsed.txt")));
		//header
		out.write("CARDSgeneShort");
		for(int i = 0; i < samples.length; i++) {
			out.write("\t" + samples[i]);
		}
		out.write("\n");
		keys = new ArrayList<String>(shortNameHits.keySet());
		Collections.sort(keys);
		for(String k : keys) {
			out.write(k);
			ArrayList<Set<String>> matches = shortNameHits.get(k);
			for(int gen = 0; gen < matches.size(); gen++) {
				if(matches.get(gen) == null || matches.get(gen).isEmpty()) {
					out.write("\tNA");
				} else {
					ArrayList<String> genes = new ArrayList<String>(matches.get(gen));
					out.write("\t" + genes.get(0));
					for(int g = 1; g < genes.size(); g++) {
						out.write(";" + genes.get(g));
					}
				}
			}
			out.write("\n");
		}
		out.close();
	}

	//if the given region overlaps an already existing region in the set, 
	//only include the region with the best hit
	//region is given as [scaffold, start, stop, bit]
	private static Set<String> checkOverlap(String region,
			Set<String> set) {
		set.add(region);
		String[] sp1 = region.split(",");
		Set<String> remove = new HashSet<String>(); //can't remove as iterating through, so remove at end
		for(String s : set) {
			String[] sp2 = s.split(",");
			//check scaffolds
			if(sp1[0].equals(sp2[0])) {
				int start1 = Integer.parseInt(sp1[1]);
				int stop1 = Integer.parseInt(sp1[2]);
				int start2 = Integer.parseInt(sp2[1]);
				int stop2 = Integer.parseInt(sp2[2]);
				if((start1 <= start2 && stop1 >= stop2) ||
						(start1 <= start2 && stop1 >= start2) ||
						(start2 <= start1 && stop2 >= start1) ||
						(start2 <= start1 && stop2 >= stop1)) {//overlap
					double bit1 = Double.parseDouble(sp1[3].replace("]", ""));
					double bit2 = Double.parseDouble(sp2[3].replace("]", ""));
					if(bit1 > bit2) {
						remove.add(s);
					} else {
						remove.add(region);
					}
				}
			}
		}
		set.removeAll(remove);
		return set;
	}

}
