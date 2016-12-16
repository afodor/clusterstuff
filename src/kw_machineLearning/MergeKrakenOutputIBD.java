/*
 * Merge the kraken outputs for China WGS, filtering human reads 
 * also generate tables split by level
 * 
 * kraken output is two columns: read classification
 */
package kw_machineLearning;

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

public class MergeKrakenOutputIBD {
	public static String DIR = "/nobackup/afodor_research/kwinglee/machineLearning/ibd/minikraken/";
	public static int NUM_SAMP;//number of samples
	public static String META = "/nobackup/afodor_research/kwinglee/machineLearning/MetAML/metaml/data/abundance_stoolsubset.txt";

	public static void main(String[] args) throws Exception {
		//get metadata
		HashMap<String, String> metaMap = new HashMap<String, String>();//map of sample ID to disease status 
		BufferedReader m = new BufferedReader(new FileReader(new File(META)));
		String[] dataset = m.readLine().split("\t");
		m.readLine();//subjectID
		String[] sampleID = m.readLine().split("\t");
		String[] disease = m.readLine().split("\t");
		m.close();
		for(int i = 0; i < dataset.length; i++) {
			if(dataset[i].equals("metahit")) {
				metaMap.put(sampleID[i], disease[i]);
			}
		}
		System.out.println("metadata " + metaMap.size());
		ArrayList<String> keys = new ArrayList<String>(metaMap.keySet());
		Collections.sort(keys);
		for(String key : keys) {
			System.out.println(key + "\t" + metaMap.get(key));
		}
		System.out.println();
		
		//get list of files to read
		ArrayList<String> tables = new ArrayList<String> ();
		String[] files = new File(DIR).list();
		for(String f : files) {
			if(f.endsWith("_mpa")) {
				tables.add(f);
			}
		}
		Collections.sort(tables);
		NUM_SAMP = tables.size();
		
		//check sequence ids
		System.out.println("sequences " + tables.size());
		HashSet<String> seqs = new HashSet<String>();
		for(int i = 0; i < tables.size(); i++) {
			String id = tables.get(i).split("_")[0].replace(".", "_");
			System.out.println(id + "\t" + metaMap.containsKey(id));
			seqs.add(id);
		}
		System.out.println();
		
		//check have all metadata
		System.out.println("missing samples");
		for(String key : keys) {
			if(!seqs.contains(key)) {
				System.out.println(key);
			}
		}
		
		
		
		/*if(tables.size() != NUM_SAMP) {
			throw new Exception("Wrong number mpa files " + tables.size());
		}

		//map of phylogeny to counts for each sample
		HashMap<String, Integer[]> baseMap = new HashMap<String, Integer[]>();

		//files set up as read phylogeny -> convert to counts
		for(int i = 0; i < tables.size(); i++) {
			BufferedReader br = new BufferedReader(new FileReader(
					new File(DIR + tables.get(i))));
			//String sample = tables.get(i).replace("minikrakenSeqs_", "").replace("_mpa", "");
			String sample = tables.get(i).replace("stdkrakenSeqs_", "").replace("_mpa", "");
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
				} else if(!sp[1].equals("root")) {
					System.err.println("human read mapped " + sample + " " + line);
					numHum++;
				} else {
					numHum++;
				}
				line = br.readLine();
			}
			br.close();
			System.out.println(sample + " " + numHum + " human reads");
		}

		////write table
		ArrayList<String> keys = new ArrayList<String>(baseMap.keySet());
		Collections.sort(keys);
		BufferedWriter out = new BufferedWriter(new FileWriter(
				//new File(DIR + "minikraken_merged.txt")));
				new File(DIR + "stdkraken_merged.txt")));
		//write header
		out.write("taxonomy");
		for(String t: tables) {
			//String sample = t.replace("minikrakenSeqs_", "").replace("_mpa", "");
			String sample = t.replace("stdkrakenSeqs_", "").replace("_mpa", "");
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
		writeSplitTable(tables, "species", split.get("species"));*/
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
				//new File(DIR + "minikraken_merged_" + level + ".txt")));
				new File(DIR + "stdkraken_merged_" + level + ".txt")));
		
		//write header
		out.write("taxa\ttaxonomy");
		for(String t: tables) {
			//String sample = t.replace("minikrakenSeqs_", "").replace("_mpa", "");
			String sample = t.replace("stdkrakenSeqs_", "").replace("_mpa", "");
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
}
