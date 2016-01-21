/**
 * From the rbh orthologs output (orthologs between two genomes),
 * for each genome generate a table of the orthologs in the other genomes
 * and a second table with the bit score of those orthologs
 */
package kw_rbh;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;

public class orthologTables implements Runnable{
	public static String DIR = "/nobackup/afodor_research/kwinglee/cre/rbh/rbhOrthologs/";
	public static String[] GROUPS = {"carolina", "resistant", "susceptible"};
	public static String ResultsFile = DIR + "orthologTables/";
	
	private String group;//group for this thread
	public orthologTables(String g) {
		group = g;
	}
	
	public static void main(String[] args) throws InterruptedException {
		new File(ResultsFile).mkdirs();
		Thread[] allThreads = new Thread[GROUPS.length];
		for(int i = 0; i < GROUPS.length; i++) {
			Runnable r = new orthologTables(GROUPS[i]);
			allThreads[i] = new Thread(r);
			allThreads[i].start();
		}
		
		//make sure all threads finish
		for(int i = 0; i < allThreads.length; i++) {
			allThreads[i].join();
		}
	}
	
	public void run() {
		File[] genomes = new File(DIR + group + "_v_" + group).listFiles();
		try {
			for(File gen : genomes) {
				String gen1 = gen.getName();
				//put all orthologs into list of maps
				List<Map<String, String>> orthName = new ArrayList<Map<String, String>>();//list of maps of gene to ortholog name
				List<Map<String, String>> bitScore = new ArrayList<Map<String, String>>();//list of maps of gene to bit score
				for(String g2 : GROUPS) {
					File[] comparisons = new File(DIR + group + "_v_" + g2 + "/" + gen1).listFiles();
					//HashMap<String, String>[] orthName = new HashMap[338];
					for(File c : comparisons) {
						String cName = c.getName().replace(".txt", "").split("_v_")[1];
						Map<String, String> names = new HashMap<String, String>();
						Map<String, String> scores = new HashMap<String, String>();
						names.put("genomeID", cName);
						scores.put("genomeID", cName);
						BufferedReader br = new BufferedReader(new FileReader(c));
						String line = br.readLine();//header
						line = br.readLine();
						while(line != null) {
							String[] sp = line.split("\t");
							names.put(sp[0], sp[1]);
							scores.put(sp[0], sp[2]);
							line = br.readLine();
						}
						br.close();
						orthName.add(names);
						bitScore.add(scores);
					}
				}
				//write results (table with ortholog name)
				printMapList(ResultsFile + "orthologNameTable_" + gen1, orthName);

				//write table with bit score
				printMapList(ResultsFile + "bitScoreTable_" + gen1, bitScore);	
			}
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	
	/*
	 * write contents of list to a file with prefix in its name
	 */
	public static void printMapList(String prefix, List<Map<String, String>> list) throws IOException {
		BufferedWriter out = new BufferedWriter(new FileWriter(new File(
				prefix + ".txt")));
		Set<String> keys = null;
		List<String> klist = null;
		for(Map<String, String> map : list) {
			Set<String> k = map.keySet(); 
			if(keys == null) {
				keys = k;
				k.remove("genomeID");
				klist = new ArrayList<String>(k);
				Collections.sort(klist);
				out.write("genomeID");
				//write header
				for(String gene : klist) {
					out.write("\t" + gene);
				}
				out.write("\n");
				System.out.println(prefix + "\t" + k.size() + "\t" + klist.size());
			} else {//check that keys are the same
				Set<String> test = k;
				test.removeAll(keys);
				if(test.size() != 0) {
					System.err.println("EXTRA IDs in:" + prefix + "/t" + map.get("genomeID"));
				}
				test = keys;
				test.removeAll(k);
				if(test.size() != 0) {
					System.err.println("EXTRA IDs in:" + prefix + "/t" + map.get("genomeID"));
				}
			}
			//write hits
			out.write(map.get("genomeID"));
			for(String gene : klist) {
				out.write("\t" + map.get(gene));
			}
			out.write("\n");
		}
		
		out.close();
	}
}
