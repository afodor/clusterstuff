/*
 * parse the results of blasting chs genomes against cards protein homolog database
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

public class ParseCardsBlastResults {
	public static int NUM_SAMPLES = 76;
	public static String DIR = "/nobackup/afodor_research/kwinglee/cre/chs_v_cards/";
	public static double PID_CUT = 80;//cutoff for percent identity
	public static double LEN_CUT = 90;//cutoff for length
	
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
		
		//get list of genes
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
		
		//for gene blasts, for each genome get set of genes that mapped to each cards gene
		HashMap<String, ArrayList<Set<String>>> hits = new HashMap<String, ArrayList<Set<String>>>();
		//HashMap<String, Integer> genomeHits = new HashMap<String, Integer>();
		HashMap<String, ArrayList<Set<String>>> shortNameHits = new HashMap<String, ArrayList<Set<String>>>();
		for(int i = 0; i < samples.length; i++) {
			String samp = samples[i];
			int numHits = 0;
			BufferedReader br = new BufferedReader(new FileReader(new File(
					DIR + "carolina_" + samp + ".genes_v_cardsProHomolog")));
			for(String line = br.readLine(); line != null; line = br.readLine()) {
				if(!line.startsWith("#")) {
					String[] sp = line.split("\t");
					String cards = sp[1];
					double pid = Double.parseDouble(sp[2]);//percent identity
					double len = 100.0 * (Double.parseDouble(sp[9]) - 
							Double.parseDouble(sp[8])) / cardsLengths.get(cards);//cards length
					if(pid > PID_CUT && len > LEN_CUT) {
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
						sets.get(i).add(sp[0]);
						
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
						sets.get(i).add(sp[0]);
					}
				}
			}
			System.out.println(samp + " " + numHits + " gene hits");
			//genomeHits.put(samp, new Integer(numHits));
			br.close();
		}

		//write results
		BufferedWriter out = new BufferedWriter(new FileWriter(new File(
				DIR + "blastCardsProHomologGeneTables_pid" + PID_CUT 
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
				if(matches.get(gen) == null) {
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
				DIR + "blastCardsProHomologGeneTables_pid" + PID_CUT 
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
				if(matches.get(gen) == null) {
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
		
		//check with scaffolds
		/*for(int i = 0; i < samples.length; i++) {
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
						numHits++;
						if(!hits.containsKey(cards)) {
							System.err.println("extra cards gene: " + cards + " in " + samp);
						}
						if(hits.get(cards) == null) {
							System.err.println("missing gene: " + cards + " in " + samp);
						}
					}
				}
			}
			System.out.println(samp + " " + numHits + " scaffold hits");
			if(numHits != genomeHits.get(samp)) {
				System.err.println("Different number hits for " + samp
						+ ": " + genomeHits.get(samp) + " genes vs "
						+ numHits + " scaffolds");
			}
			br.close();
		}*/
	}

}
