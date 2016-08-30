/*
 * for each of the genes identified as a beta lactamase, identify which class
 * is the best hit
 * 8/5/16
 */

package kw_cre;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;

public class GetBestBlastHits_nonBetaLactamases {
	public static String DIR = "/nobackup/afodor_research/kwinglee/cre/chs_v_cards/";

	public static void main(String[] args) throws IOException {
		String[] files = new File(DIR).list();
		HashMap<String, String> bitScore = new HashMap<String, String>();//gene-> best hit based on bit score (hit saved as bit score; percent identity; gene)
		HashMap<String, String> pID = new HashMap<String, String>();//gene-> best hit based on percent identity (hit saved as bit score; percent identity; gene)
		for(String f : files) {
			if(f.endsWith("genes_v_cardsProHomolog")) {
				BufferedReader blast = new BufferedReader(new FileReader(new File(
						DIR + f)));
				for(String line = blast.readLine(); line != null; line = blast.readLine()) {
					if(!line.startsWith("#")) {
						String[] sp = line.split("\t");
						String gene = sp[0];
						String cards = sp[1];
						double pid = Double.parseDouble(sp[2]);//percent identity
						double bit = Double.parseDouble(sp[11]);//bit score
						String len = sp[3];
						String value = bit + ";" + pid + ";" + len + ";" + cards;
						if(bitScore.containsKey(gene)) {
							//see if best bit score
							String bestBit = bitScore.get(gene);
							double b = Double.parseDouble(bestBit.split(";")[0]);
							if(bit > b) {
								bitScore.put(gene, value);
							} else if (bit == b) {
								double p = Double.parseDouble(bestBit.split(";")[1]);
								if(pid > p) {
									bitScore.put(gene, value);
								} else if(pid == p) {
									bitScore.put(gene, bestBit + "," + cards);
								}
							}
							
							//see if best percent identity
							String bestPid = pID.get(gene);
							double p = Double.parseDouble(bestPid.split(";")[1]);
							if(pid > p) {
								pID.put(gene, value);
							} else if(pid == p) {
								b = Double.parseDouble(bestPid.split(";")[0]);
								if(bit > b) {
									pID.put(gene, value);
								} else if(bit == b) {
									pID.put(gene, bestPid + "," + cards);
								}
							}
						} else {
							bitScore.put(gene, value);
							pID.put(gene, value);
						}
					}
				}
				blast.close();
			}
		}
		
		////make map of collapsed CARD gene names to hits
		HashMap<String, ArrayList<String>> hits = new HashMap<String, ArrayList<String>>();//map of collapsed card name to all genes in that mapped
		BufferedReader br = new BufferedReader(new FileReader(new File(
				DIR + "blastCardsProHomologGeneTables_pid80.0_len90.0_collapsed.txt")));
		br.readLine();//header
		for(String line = br.readLine(); line != null; line = br.readLine()) {
			String[] sp = line.split("\t");
			String key = sp[0];
			ArrayList<String> list = new ArrayList<String>(Arrays.asList(sp));
			list.remove(key);
			hits.put(key, list);
		}
		br.close();
		
		////look specifically at certain genes
		String[] list = new String[]{"acrA", "acrB", "acrD", "APH(3'')", "APH(6)", 
				"emrB", "emrD", "mdtB", "mdtC", "mdtD", "sul1"};
		for(String card : list) {
			BufferedWriter out = new BufferedWriter(new FileWriter(new File(
					DIR + "bestBlastHit_nonBetaLactamase_" + card + ".txt")));
			out.write("gene\tbestCardsHit\tbitScore\tpercentID\tlength\tbestMatchForThisCard\n");
			ArrayList<String> matches = hits.get(card);
			for(String m : matches) {
				String[] sp = m.split(";");
				for(String s : sp) {
					if(!s.equals("NA")) {
						String[] bit = bitScore.get(s).split(";");
						
						out.write(s + "\t" + bit[3] + "\t" + bit[0] + "\t" 
								+ bit[1] + "\t" + bit[2] + "\t" +
								bit[3].contains(card) + "\n");
						
					}
				}
			}
			out.close();
		}

	}
}
