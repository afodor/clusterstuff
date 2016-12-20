/*
 * Merge the kraken outputs for IBD 
 * also generate tables split by level
 * 
 * kraken output is two columns: read classification
 */
package kw_jobinMiRNA;

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
import java.util.concurrent.Semaphore;

public class MergeKrakenOutput {
	public static String DIR = "/nobackup/afodor_research/kwinglee/jobin/microRNA/";

	public static void main(String[] args) throws InterruptedException {
		//set up multithreaded
		String[] names = new String[]{"miniKraken", "stdKraken"};
		Semaphore s = new Semaphore(names.length);

		for(String n : names) {
			s.acquire();
			MergeOutput f = new MergeOutput(n, s);
			new Thread(f).start();
		}

		for(int i = 0; i < names.length; i++) {
			s.acquire();
		}
	}
}

class MergeOutput implements Runnable {
	private String name;
	private final Semaphore semaphore;
	private final static int NUM_SAMP = 10;


	public MergeOutput(String name, Semaphore s) {
		this.name = name;
		this.semaphore = s;
	}


	@Override
	public void run() {
		try {
			//get list of files to read
			ArrayList<String> tables = new ArrayList<String> ();
			String[] files = new File(MergeKrakenOutput.DIR + name).list();
			for(String f : files) {
				if(f.endsWith("_mpa")) {
					tables.add(f);
				}
			}
			Collections.sort(tables);

			//map of phylogeny to counts for each sample
			HashMap<String, Integer[]> baseMap = new HashMap<String, Integer[]>();

			//files set up as read phylogeny -> convert to counts
			for(int tab = 0; tab < tables.size(); tab++) {
				BufferedReader br = new BufferedReader(new FileReader(
						new File(MergeKrakenOutput.DIR +
								name + File.separator + tables.get(tab))));
				String line = br.readLine();
				while(line != null) {
					String[] sp = line.split("\t");
					String taxa = sp[1];
					if(!baseMap.containsKey(taxa)) {//taxa not previously seen; add array of zeros
						Integer[] counts = new Integer[NUM_SAMP];
						Arrays.fill(counts, 0);
						baseMap.put(taxa, counts);
					}
					Integer[] counts = baseMap.get(taxa);
					counts[tab]++;
					line = br.readLine();
				}
				br.close();
			}

			////write table
			ArrayList<String> keys = new ArrayList<String>(baseMap.keySet());
			Collections.sort(keys);
			BufferedWriter out = new BufferedWriter(new FileWriter(
					new File(MergeKrakenOutput.DIR + name
							+ File.separator + "miRNA_" + name.toLowerCase() + ".txt")));
			//write header
			out.write("taxonomy");
			for(String t: tables) {
				out.write("\t" + t.replace("_mpa", ""));
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
		} catch (Exception e) {
			System.err.println(e.getMessage());
		} finally {
			semaphore.release();
		}
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
	private void writeSplitTable(ArrayList<String> ids, 
			String level, 
			HashMap<String, Integer[]> map) throws IOException {
		BufferedWriter out = new BufferedWriter(new FileWriter(
				new File(MergeKrakenOutput.DIR + name
						+ File.separator + "miRNA_" + name.toLowerCase() + 
						"_" + level + ".txt")));

		//write header
		out.write("taxa\ttaxonomy");
		for(String t: ids) {
			out.write("\t" + t.replace("_mpa", ""));
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
