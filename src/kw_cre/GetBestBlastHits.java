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
import java.util.HashMap;

public class GetBestBlastHits {
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
		
		////look specifically at the beta lactamase genes
		String bldir = DIR + "betaLactamaseAlignments/";
		String[] blList = new String[]{"KPC", "LEN", "OKP", "OXA", "SHV", "TEM"};
		ArrayList<String> seenGenes = new ArrayList<String>();
		
		for(String bl : blList) {
			//get list of genes for that betalactamase
			ArrayList<String> genes = new ArrayList<String>();
			BufferedReader fasta = new BufferedReader(new FileReader(new File(
					bldir + bl + ".fasta")));
			for(String line = fasta.readLine(); line != null; line = fasta.readLine()) {
				if(line.startsWith(">")) {
					String gene = line.replace(">", ""); 
					genes.add(gene);
					if(seenGenes.contains(gene)) {
						System.out.println(bl + "\t" + gene);
					} else {
						seenGenes.add(gene);
					}
				}
			}
			fasta.close();
			
			//print best hit (checking that the best percent identity and bit score are the same)
			BufferedWriter out = new BufferedWriter(new FileWriter(new File(
					DIR + "bestBlastHit_" + bl + ".txt")));
			out.write("gene\tbestCardsHit\tbitScore\tpercentID\tlength\n");
			BufferedWriter alt = new BufferedWriter(new FileWriter(new File(
					DIR + "bestBlastHit_" + bl + "_diffBestPID.txt")));
			alt.write("gene\tbestCardsHit\tbitScore\tpercentID\tlength\n");
			for(String g : genes) {
				if(bitScore.containsKey(g)) {
					String value = bitScore.get(g);
					String[] sp = value.split(";");
					out.write(g + "\t" + sp[3] + "\t" + sp[0] + "\t" + sp[1] + "\t" + sp[2] + "\n");
					if(!value.equals(pID.get(g))) {
						value = pID.get(g);
						sp = value.split(";");
						alt.write(g + "\t" + sp[3] + "\t" + sp[0] + "\t" + sp[1] + "\t" + sp[2] + "\n");
						//System.out.println(g + "\t" + pID.get(g));
					}
				} else {
					System.err.println(g + " is missing from files");
				}
			}
			out.close();
			alt.close();
		}
	}
}
