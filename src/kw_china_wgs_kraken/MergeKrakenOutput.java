/*
 * Merge the kraken outputs for China WGS, filtering human reads 
 * also generate tables split by level
 * 
 * kraken output is two columns: read classification
 */
package kw_china_wgs_kraken;

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
import java.util.Set;

public class MergeKrakenOutput {
	public static String DIR = "/nobackup/afodor_research/kwinglee/china/wgs/minikrakenResults/";
	public static String HGDIR= "/nobackup/afodor_research/kwinglee/china/wgs/alignToHG38/";
	public static int NUM_SAMP = 40;//number of samples

	public static void main(String[] args) throws Exception {
		//get list of files to read
		ArrayList<String> tables = new ArrayList<String> ();
		String[] files = new File(DIR).list();
		for(String f : files) {
			if(f.endsWith("_mpa")) {
				tables.add(f);
			}
		}
		Collections.sort(tables);
		if(tables.size() != NUM_SAMP) {
			throw new Exception("Wrong number mpa files " + tables.size());
		}

		//map of phylogeny to counts for each sample
		HashMap<String, Integer[]> baseMap = new HashMap<String, Integer[]>();

		//files set up as read phylogeny -> convert to counts
		for(int i = 0; i < tables.size(); i++) {
			BufferedReader br = new BufferedReader(new FileReader(
					new File(DIR + tables.get(i))));
			String sample = tables.get(i).replace("minikrakenSeqs_", "").replace("_mpa", "");
			Set<String> human = getHumanReads(sample);
			int numHum = 0;//number human reads seen
			String line = br.readLine();
			while(line != null) {
				String[] sp = line.split("\t");
				String taxa = sp[1];
				if(!human.contains(sp[0])) {
					if(!baseMap.containsKey(taxa)) {//taxa not previously seen; add array of zeros
						Integer[] counts = new Integer[tables.size()];
						Arrays.fill(counts, 0);
						baseMap.put(taxa, counts);
					}
					Integer[] counts = baseMap.get(taxa);
					counts[i]++;
					line = br.readLine();
				} else if(!sp[1].equals("root")) {
					br.close();
					throw new Exception("human read mapped " + sample + " " + sp[0]);
				} else {
					numHum++;
				}
			}
			br.close();
			System.out.println(sample + " " + numHum + "human reads");
		}

		////write table
		ArrayList<String> keys = new ArrayList<String>(baseMap.keySet());
		Collections.sort(keys);
		BufferedWriter out = new BufferedWriter(new FileWriter(
				new File(DIR + "minikraken_merged.txt")));
		//write header
		out.write("taxonomy");
		for(String t: tables) {
			String sample = t.replace("minikrakenSeqs_", "").replace("_mpa", "");
			out.write("\t" + sample);
		}
		out.write("\n");
		//write counts
		for(String k : keys) {
			out.write(k);
			Integer[] counts = baseMap.get(k);
			for(int i = 0; i < counts.length; i++) {
				out.write("\t" + counts[i]);
			}
			out.write("\n");
		}
		out.close();

		////split by level
		//map of level -> (map of taxa -> counts)
		HashMap<String, HashMap<String, Integer[]>> split = 
				new HashMap<String, HashMap<String, Integer[]>>();
		split.put("domain", new HashMap<String, Integer[]>());
		split.put("phylum", new HashMap<String, Integer[]>());
		split.put("class", new HashMap<String, Integer[]>());
		split.put("order", new HashMap<String, Integer[]>());
		split.put("family", new HashMap<String, Integer[]>());
		split.put("genus", new HashMap<String, Integer[]>());
		split.put("species", new HashMap<String, Integer[]>());

		//for each key, split by level and add counts
		for(String k : keys) {
			if(!k.equals("root")) {
				Integer[] counts = baseMap.get(k);
				String[] sp = k.split("\\|");//Pattern.quote("|")
				String name = sp[0];
				addCounts(split.get("domain"), name, counts);

				for(int i = 1; i < sp.length; i++) {
					name += "|" + sp[i];
					if(sp[i].startsWith("p__")) {
						addCounts(split.get("phylum"), name, counts);
					} else if(sp[i].startsWith("c__")) {
						addCounts(split.get("class"), name, counts);
					} else if(sp[i].startsWith("o__")) {
						addCounts(split.get("order"), name, counts);
					} else if(sp[i].startsWith("f__")) {
						addCounts(split.get("family"), name, counts);
					} else if(sp[i].startsWith("g__")) {
						addCounts(split.get("genus"), name, counts);
					} else if(sp[i].startsWith("s__")) {
						addCounts(split.get("species"), name, counts);
					} else {
						System.out.println("unknown phylogeny: " + sp[i] + " " + k);
					}
				}
			}
		}

		//write tables
		writeSplitTable(tables, "domain", split.get("domain"));
		writeSplitTable(tables, "phylum", split.get("phylum"));
		writeSplitTable(tables, "class", split.get("class"));
		writeSplitTable(tables, "order", split.get("order"));
		writeSplitTable(tables, "family", split.get("family"));
		writeSplitTable(tables, "genus", split.get("genus"));
		writeSplitTable(tables, "species", split.get("species"));
	}

	//function that adds the given counts to the appropriate key in the given map
	public static void addCounts(HashMap<String, Integer[]> map,
			String key, Integer[] newCounts) {
		if(!map.containsKey(key)) {
			Integer[] counts = new Integer[NUM_SAMP];
			Arrays.fill(counts, 0);
			map.put(key, counts);
		}
		Integer[] counts = map.get(key);
		for(int i = 0; i < counts.length; i++) {
			counts[i] += newCounts[i];
		}
	}

	//writes the table for the given level containing the given counts
	public static void writeSplitTable(ArrayList<String> tables, 
			String level, 
			HashMap<String, Integer[]> map) throws IOException {
		BufferedWriter out = new BufferedWriter(new FileWriter(
				new File(DIR + "minikraken_merged_" + level + ".txt")));

		//write header
		out.write("taxa\ttaxonomy");
		for(String t: tables) {
			String sample = t.replace("minikrakenSeqs_", "").replace("_mpa", "");
			out.write("\t" + sample);
		}
		out.write("\n");

		//write counts
		ArrayList<String> keys = new ArrayList<String>(map.keySet());
		Collections.sort(keys);
		System.out.println(level + " " + keys.size());
		for(String k : keys) {
			String[] ksplit = k.split("\\|");
			out.write(ksplit[ksplit.length-1].replaceFirst(".__", "") +
					"\t" + k);
			Integer[] counts = map.get(k);
			for(int i = 0; i < counts.length; i++) {
				out.write("\t" + counts[i]);
			}
			out.write("\n");
		}

		out.close();
	}
	
	//for the given genome, return the set of reads that mapped to the human genome
	public static Set<String> getHumanReads(String sample) throws IOException {
		Set<String> reads = new HashSet<String>();
		BufferedReader map = new BufferedReader(new FileReader(new File(
				HGDIR + sample + "_1.hg38.mapped.sam")));
		String line = map.readLine();
		while(line != null) {
			String[] sp = line.split("\t");
			reads.add(sp[0] + "/1");//kegg results and fasta have extra /1
			line = map.readLine();
		}
		map.close();
		return(reads);
	}
}
