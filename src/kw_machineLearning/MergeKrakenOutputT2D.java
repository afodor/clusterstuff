/*
 * Merge the kraken outputs for IBD 
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

public class MergeKrakenOutputT2D {
	public static String DIR = "/nobackup/afodor_research/kwinglee/machineLearning/t2d/minikraken/";
	public static int NUM_SAMP;//number of samples
	public static String META = "/nobackup/afodor_research/kwinglee/machineLearning/MetAML/metaml/data/abundance_stoolsubset.txt";
	public static String BASEOUT = "t2d_minikraken_merged";//prefix of output files
	public static String SRADIR = "/nobackup/afodor_research/kwinglee/machineLearning/t2d/";

	public static void main(String[] args) throws Exception {
		//get metadata from Segata paper
		/*HashMap<String, String> metaMap = new HashMap<String, String>();//map of sample ID to disease status 
		BufferedReader m = new BufferedReader(new FileReader(new File(META)));
		String[] dataset = m.readLine().split("\t");
		String[] sampleID = m.readLine().split("\t");
		m.readLine();//subjectID
		m.readLine();//bodysite
		String[] disease = m.readLine().split("\t");
		m.close();
		for(int i = 0; i < dataset.length; i++) {
			if(dataset[i].equals("t2dmeta_short") || dataset[i].equals("t2dmeta_long")) {
				metaMap.put(sampleID[i], disease[i]);
			}
		}
		System.out.println("metadata " + metaMap.size());
		ArrayList<String> keys = new ArrayList<String>(metaMap.keySet());
		Collections.sort(keys);*/
		/*for(String key : keys) {
			System.out.println(key + "\t" + metaMap.get(key));
		}
		System.out.println();*/
		
		//get insert size from T2D paper
		HashMap<String, String> insertMap = new HashMap<String, String>();//map of sample ID to insert size
		BufferedReader s2 = new BufferedReader(new FileReader(new File(
				SRADIR + "SuppTableS2.txt")));
		s2.readLine();//header
		s2.readLine();//header
		for(String line = s2.readLine(); line != null; line = s2.readLine()) {
			String[] sp = line.split("\t");
			if(sp.length > 2) {
				insertMap.put(sp[1].replace(" ", "-"), sp[2]);	
			}
		}
		s2.close();

		//get metadata from T2D paper
		HashMap<String, String> gaToID = new HashMap<String, String>();//map of gender/age to sample ID
		HashMap<String, String> idToGroup = new HashMap<String, String>();//map of sample ID to group
		BufferedReader pprTab = new BufferedReader(new FileReader(new File(
				SRADIR + "SuppTableS1.txt")));
		pprTab.readLine();//header
		for(String line = pprTab.readLine(); line != null; line = pprTab.readLine()) {
			String[] sp = line.split("\t");
			String id = sp[1].replace(" ", "-");
			String ga = sp[2] + sp[3] + insertMap.get(id);//gender+age+insert size
			ga = ga.replaceAll("\\s","");
			String group = sp[7];
			if(gaToID.containsKey(ga)) {
				gaToID.put(ga, gaToID.get(ga) + ";" + id);
			} else {
				gaToID.put(ga, id);				
			}

			if(group.equals("N")) {
				group = "n";
			} else {
				group = "t2d";
			}
			idToGroup.put(id, group);
		}
		pprTab.close();
		/*System.out.println("ga " + gaToID.size() + " id " + idToGroup.size() + 
				" insert " + insertMap.size());*/

		//check Segata and paper tables are giving similar results
		/*for(String k : keys) {
			if(idToGroup.containsKey(k)) {
				if(!idToGroup.get(k).equals(metaMap.get(k))) {
					System.out.println("Different class: " + k + " " 
							+ idToGroup.get(k) + " " + metaMap.get(k));
				}
			} else {
				System.out.println("Missing sample: " + k + " " + metaMap.get(k));
			}
		}*/

		//get list of files to read
		ArrayList<String> tables = new ArrayList<String> ();
		String[] files = new File(DIR).list();
		for(String f : files) {
			if(f.endsWith("_mpa")) {
				tables.add(f);
			}
		}
		//Collections.sort(tables);

		//get conversion from SRA to sample ID
		//HashMap<String, String> srrToGA = new HashMap<String, String>();//srr to gender + age
		//HashMap<String, String> srrToID = new HashMap<String, String>();//srr to sequence sample name
		HashMap<String, String> metaMap = new HashMap<String, String>();//srr to class
		String[] sraFiles = new String[] {"SraRunTable.txt", "SraRunTable2.txt"};
		HashSet<String> matchedKeys = new HashSet<String>();//set of sample IDs seen
		HashSet<String> unmatchedSplit = new HashSet<String>();//set of sample IDs not matched due to multiple same key
		//int numMissingGaToIDkey = 0;
		//int numMult = 0;
		BufferedWriter out = new BufferedWriter(new FileWriter(new File(
				DIR + BASEOUT + "_metadataFull.txt")));
		out.write("SRR_ID\tBGI_ID\tpaper_ID\tgender\tage\tinsert_size\tgroup\n");
		for(String s : sraFiles) {
			BufferedReader br = new BufferedReader(new FileReader(new File(SRADIR + s)));
			String[] head = br.readLine().split("\t");//header
			int nameCol = 0;
			int srrCol = 0;
			int genderCol = 0;
			int ageCol = 0;
			int insertCol = 0;
			for(int i = 0; i < head.length; i++) {
				if(head[i].equals("Sample_Name_s")) {
					nameCol = i;
				} else if (head[i].equals("Run_s")) {
					srrCol = i;
				} else if(head[i].equals("GENDER_s")) {
					genderCol = i;
				} else if(head[i].equals("AGE_s")) {
					ageCol = i;
				} else if(head[i].equals("InsertSize_l")) {
					insertCol = i;
				}
			}
			/*System.out.println("Columns: " + nameCol + " " + " " + srrCol + " " + 
					genderCol + " " + ageCol);*/
			for(String line = br.readLine(); line != null; line = br.readLine()) {
				String[] sp = line.split("\t");
				if(sp.length > 1) {
					String srr = sp[srrCol];
					String ga = sp[genderCol] + sp[ageCol] + sp[insertCol];
					String seqID = sp[nameCol].replace("bgi-", "").replace(" ", "-");
					String pprID = gaToID.get(ga);
					matchedKeys.add(pprID);
					/*if(!gaToID.containsKey(ga)) {
						System.out.println("Missing gaToID key " + srr + " " + seqID
								+ " " + pprID + " " + ga);
						numMissingGaToIDkey++;
					} else {*/
					if(gaToID.containsKey(ga)) {
						if(pprID.contains(";")) {
							if(pprID.contains(seqID)) {
								String split = pprID.replace(seqID + ";", "").replace(";" + seqID, "");
								//System.out.println(seqID + " " + pprID + " " + split);
								gaToID.put(ga, split);
								pprID = seqID;
							} else {
								/*System.out.println("To split " + ga + " " +
										seqID + " " + pprID);	*/
								unmatchedSplit.add(srr + " " + ga + " " + seqID);
								//numMult++;
							}
						} /*else {
							System.out.println("Test: " + srr + " " + pprID + " " + ga);
						}*/
						if(idToGroup.containsKey(pprID)) {
							metaMap.put(srr, idToGroup.get(pprID));
						} /*else if(!pprID.contains(";")){
							System.out.println("Missing idToGroup key " + srr + " " + seqID
									+ " " + pprID + " " + ga);
						}*/
					} else if(seqID.startsWith("T2D") || seqID.startsWith("ED") ||
							seqID.startsWith("HT")) {
						metaMap.put(srr, "T2D");
					} else if(seqID.startsWith("SZEY") || seqID.startsWith("N")) {
						metaMap.put(srr, "N");
					} else {
						System.out.println(seqID);
					}
					out.write(srr + "\t" + seqID + "\t" + pprID + "\t" + sp[genderCol] + "\t"
							+ sp[ageCol] + "\t" + sp[insertCol] + "\t" + idToGroup.get(pprID) 
							+ metaMap.get(srr) + "\n");
				} /*else {
					System.out.println(line);
				}*/
			}
			br.close();
		}
		out.close();
		/*System.out.println("Number missing gaToID key: " + numMissingGaToIDkey);
		System.out.println("Number multiple same gender/age/insert key first pass: " + numMult
				+ " " + unmatchedSplit.size());*/
		
		//see if any singletons after splitting samples with same gender/age/insert size
		//numMult = 0;
		for(String s : unmatchedSplit) {
			String[] sp = s.split(" ");
			String srr = sp[0];
			String ga = sp[1];
			String seqID = sp[2];
			String pprID = gaToID.get(ga);
			matchedKeys.add(pprID);
			if(pprID.contains(";")) {
				if(pprID.contains(seqID)) {
					gaToID.put(ga, 
							pprID.replace(seqID + ";", "").replace(";" + seqID, ""));
					pprID = seqID;
				} /*else {
					System.out.println("To split " + ga + " " +
							seqID + " " + pprID);
					numMult++;
				}*/
			}
			if(idToGroup.containsKey(pprID)) {
				metaMap.put(srr, idToGroup.get(pprID));
			} 
		}
		//System.out.println("Number multiple same gender/age/insert key second pass: " + numMult);
		ArrayList<String> keys = new ArrayList<String>(metaMap.keySet());
		Collections.sort(keys);

		//check sequence ids
		System.out.println("SRA " + metaMap.size());
		//System.out.println("sequences " + tables.size());
		HashSet<String> seqs = new HashSet<String>();
		int numMissing = 0;
		for(int i = 0; i < tables.size(); i++) {
			String id = tables.get(i).split("_")[0];
			if(!metaMap.containsKey(id)) {
				//System.out.println("Missing table " + id);
				numMissing++;
			} /*else {
				System.out.println(id + "\t" + metaMap.containsKey(id) + "\t" + metaMap.get(id));
			}*/
			seqs.add(id);
		}
		//System.out.println();
		System.out.println("Number missing tables: " + numMissing);

		//check have all metadata
		//System.out.println("missing samples");
		//System.out.println("Number paper samples " + insertMap.size());
		/*numMissing = 0;
		for(String key : insertMap.keySet()) {
			if(!matchedKeys.contains(key)) {
				//System.out.println(key);
				numMissing++;
			}
		}
		System.out.println("Number missing paper samples " + numMissing);*/

		//add missing samples and make map of sample id to sequences
		/*HashMap<String, ArrayList<String>> sequences = new HashMap<String, ArrayList<String>>();//map of id to the sequences associated 
		for(int i = 0; i < tables.size(); i++) {
			String id = tables.get(i).split("_")[0];
			if(!sraMap.containsKey(id) || !metaMap.containsKey(sraMap.get(sra))) {
				System.out.println("Extra key: " + id);
			}
			if(sequences.containsKey(id)) {
				sequences.get(id).add(tables.get(i));
			} else {
				ArrayList<String> list = new ArrayList<String>();
				list.add(tables.get(i));
				sequences.put(id, list);
			}
		}
		ArrayList<String> seqIDs = new ArrayList<String>(sequences.keySet());
		Collections.sort(seqIDs);
		NUM_SAMP = seqIDs.size();

		//map of phylogeny to counts for each sample
		HashMap<String, Integer[]> baseMap = new HashMap<String, Integer[]>();

		//files set up as read phylogeny -> convert to counts
		for(int key = 0; key < seqIDs.size(); key++) {
			ArrayList<String> tabs = sequences.get(seqIDs.get(key));
			for(int tab = 0; tab < tabs.size(); tab++) {
				BufferedReader br = new BufferedReader(new FileReader(
						new File(DIR + tabs.get(tab))));
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
					counts[key]++;
					line = br.readLine();
				}
				br.close();
			}
		}

		////write table
		ArrayList<String> keys = new ArrayList<String>(baseMap.keySet());
		Collections.sort(keys);
		BufferedWriter out = new BufferedWriter(new FileWriter(
				new File(DIR + BASEOUT + ".txt")));
		//write header
		out.write("taxonomy");
		for(String t: seqIDs) {
			out.write("\t" + t + "." + sraMap.get(t));//id is SRA.SampleName
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
		writeSplitTable(seqIDs, "domain", split.get("domain"), sraMap);
		writeSplitTable(seqIDs, "phylum", split.get("phylum"), sraMap);
		writeSplitTable(seqIDs, "class", split.get("class"), sraMap);
		writeSplitTable(seqIDs, "order", split.get("order"), sraMap);
		writeSplitTable(seqIDs, "family", split.get("family"), sraMap);
		writeSplitTable(seqIDs, "genus", split.get("genus"), sraMap);
		writeSplitTable(seqIDs, "species", split.get("species"), sraMap);

		//write metadata table
		out = new BufferedWriter(new FileWriter(new File(
				DIR + BASEOUT + "_metadata.txt")));
		out.write("sampleID\tdisease\n");
		for(String s : seqIDs) {
			out.write(s + "." + sraMap.get(s) + "\t" 
					+ metaMap.get(sraMap.get(s)) + "\n");
		}
		out.close();*/
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
	public static void writeSplitTable(ArrayList<String> ids, 
			String level, 
			HashMap<String, Integer[]> map,
			HashMap<String, String> sraMap) throws IOException {
		BufferedWriter out = new BufferedWriter(new FileWriter(
				new File(DIR + BASEOUT + "_" + level + ".txt")));

		//write header
		out.write("taxa\ttaxonomy");
		for(String t: ids) {
			out.write("\t" + t + "." + sraMap.get(t));//id is SRA.SampleName
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
