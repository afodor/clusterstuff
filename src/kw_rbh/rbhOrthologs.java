/**
 * From BLAST results, identify RBH orthologs:
 * reciprocal best blast hit are those where for sequence sa in genome A and 
 * sequence sb in genome B the best match to sa is sb and the best match to sb is sa
 * 
 * output is table for each genome comparison where for each gene in the first genome the 
 * ortholog in the other genome is given, with the bit scores for each
 */
package kw_rbh;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

public class rbhOrthologs implements Runnable {
	public static String DIR = "/nobackup/afodor_research/kwinglee/cre/rbh/";
		
	public static void main(String[] args) throws InterruptedException {
		Thread[] allThreads = new Thread[9];
		Runnable r = new rbhOrthologs("carolina_v_carolina", "carolina_v_carolina");
		allThreads[0] = new Thread(r);
		r = new rbhOrthologs("carolina_v_resistant", "resistant_v_carolina");
		allThreads[1] = new Thread(r);
		r = new rbhOrthologs("carolina_v_susceptible", "susceptible_v_carolina");
		allThreads[2] = new Thread(r);
		r = new rbhOrthologs("resistant_v_resistant", "resistant_v_resistant");
		allThreads[3] = new Thread(r);
		r = new rbhOrthologs("resistant_v_carolina", "carolina_v_resistant");
		allThreads[4] = new Thread(r);
		r = new rbhOrthologs("resistant_v_susceptible", "susceptible_v_resistant");
		allThreads[5] = new Thread(r);
		r = new rbhOrthologs("susceptible_v_susceptible", "susceptible_v_susceptible");
		allThreads[6] = new Thread(r);
		r = new rbhOrthologs("susceptible_v_carolina", "carolina_v_susceptible");
		allThreads[7] = new Thread(r);
		r = new rbhOrthologs("susceptible_v_resistant", "resistant_v_susceptible");
		allThreads[8] = new Thread(r);
		
		for(int i = 0; i < allThreads.length; i++) {
			allThreads[i].start();
		}
		
		//make sure all threads finish
		for(int i = 0; i < allThreads.length; i++) {
			allThreads[i].join();
		}
	}
	
	/**
	 * 
	 * @param folder1	folder with set of genomes
	 * @param folder2	folder with reciprocal set of genomes
	 * @param outFolder	folder to write results to
	 */
	private String folder1;
	private String folder2;
	public rbhOrthologs(String f1, String f2) {
		folder1 = f1;
		folder2 = f2;
		new File(DIR + "rbhOrthologs/" + folder1).mkdirs();
	}
	
	public void run() {
		File f1 = new File(DIR + "blastResults/" + folder1);
		File[] list1 = f1.listFiles();
		for(File f : list1) {
			File[] blasts = f.listFiles();
			new File(DIR + "rbhOrthologs/" + folder1 + "/" + f.getName()).mkdirs();
			for(File b : blasts) {
				try {
					String[] name = b.getName().replace(".txt", "").split("_v_");
					//if(!name[0].equals(name[1])) {//only analyze if not comparing to self
						HashMap<String, String[]> map1 = getHits(b);
						HashMap<String, String[]> map2 = getHits(new File(
								DIR + "blastResults/" + folder2 + "/" + name[1] + "/" + name[1] + "_v_" + name[0] + ".txt"));
						List<String> genes = getGenes(name[0]);
						BufferedWriter out = new BufferedWriter(new FileWriter(new File(
								DIR + "rbhOrthologs/" + folder1 + "/" + f.getName() + "/rbhResults_" + b.getName())));
						out.write("Gene\tortholog\tBitScore1\tBitScore2\n");
						for(String g : genes) {
							if(map1.containsKey(g)) {
								String[] match = map1.get(g);
								String g2 = match[0];
								if(map2.containsKey(g2)) {
									String[] match2 = map2.get(g2);
									if(match2[0].equals(g)) {//gene is reciprocal
										out.write(g + "\t" + g2 + "\t" + match[1] + "\t" + match2[1] + "\n");
									} else {
										out.write(g + "NA\t0\t0\n");
									}
								} else {
									out.write(g + "NA\t0\t0\n");
								}
							} else { //gene is not hit it both genomes
								out.write(g + "NA\t0\t0\n");
							}
						}
						out.close();
					//}
				} catch (IOException e) {
					System.err.println("ERROR IN: " + folder1);
					e.printStackTrace();
					break;
				}
			}
		}
	}
	
	//for the given genome, return list of all genes
	public static List<String> getGenes(String genome) throws IOException {
		String[] name = genome.split("_");
		List<String> genes = new ArrayList<String>();
		BufferedReader br = new BufferedReader(new FileReader(new File(
				DIR + name[0] + "/" + genome + "_allGenes.fasta")));
		String line = br.readLine();
		while(line != null) {
			if(line.startsWith(">")) {
				genes.add(line.replace(">", ""));
			}
			line = br.readLine();
		}
		br.close();
		return(genes);
	}
	
	//for given blast results file, return hash of gene name to top hit and its bit score
	public static HashMap<String, String[]> getHits(File file) throws IOException {
		HashMap<String, String[]> map = new HashMap<String, String[]>();
		BufferedReader br = new BufferedReader(new FileReader(file));
		String line = br.readLine();
		while(line != null) {
			if(!line.startsWith("#")) {
				String[] sp = line.split("\t");
				String key = sp[0].trim();
				String match = sp[1].trim();
				String score = sp[11].trim();
				if(!map.containsKey(key)) {
					map.put(key, new String[] {match, score});
				} else {//key is already in there
					String prevScore = map.get(key)[1];
					if(Double.parseDouble(score) > Double.parseDouble(prevScore)) {//replace with better hit
						map.put(key, new String[] {match, score});
					}
				}
			}
			line = br.readLine();
		}
		br.close();
		return(map);
	}
}
