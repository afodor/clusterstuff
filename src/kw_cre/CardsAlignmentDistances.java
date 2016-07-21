/*
 * calculates the distances between genes from the alignment resulting from
 * AlignSelectCardsHits.java
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
import java.util.List;

public class CardsAlignmentDistances {
	public static String DIR = "/nobackup/afodor_research/kwinglee/cre/chs_v_cards/betaLactamaseAlignments/";
	
	public static void main(String[] args) throws Exception {
		String[] betalacs = new String[]{"KPC", "LEN", "OKP", "OXA", "SHV", "TEM"};//genes to analyze
		//String[] betalacs = new String[]{"test"};
		
		for(String bl : betalacs) {
			HashMap<String, String> seqs = new HashMap<String, String>();//map of gene to aligned sequence
			//get sequences
			BufferedReader fasta = new BufferedReader(new FileReader(new File(
					DIR + bl + ".mafftAlign.fasta")));
			String header = "";
			String align = "";
			for(String line = fasta.readLine(); line!= null; line = fasta.readLine()) {
				if(line.startsWith(">")) {
					if(header.length() > 1) {
						seqs.put(header, align);
					}
					header = line.replace(">", "");
					align = "";
				} else {
					align += line;
				}
			}
			seqs.put(header, align);
			fasta.close();
			List<String> genes = new ArrayList<String>();
			genes.addAll(seqs.keySet());
			Collections.sort(genes);
			
			//set up distance matrix
			int[][] dist = new int[genes.size()][];
			for(int i = 0; i < dist.length; i++) {
				dist[i] = new int[genes.size()];
			}
			
			//calculate distances
			boolean nonzero = false;
			for(int i = 0; i < genes.size(); i++) {
				String gene1 = genes.get(i);
				String seq1 = seqs.get(gene1);
				for(int j = 0; j < genes.size(); j++) {
					String gene2 = genes.get(j);
					String seq2 = seqs.get(gene2);
					if(seq1.length() != seq2.length()) {
						throw new Exception("Different lengths " + gene1 
								+ " " + gene2);
					} else {
						int diff = 0;
						boolean gap = false;
						for(int k = 0; k < seq1.length(); k++) {
							if(seq1.charAt(k) != (seq2.charAt(k)) && !gap) {//sequence is different and not currently part of gap
								diff++;
								/*if(bl == "KPC") {
									System.out.println(gene1 + "\t" + gene2 + "\t" + j);
								}*/
								gap = seq1.charAt(k) == '-' || seq2.charAt(k) == '-';
							} else if(gap) {//possible end of gap
								gap = seq1.charAt(k) == '-' || seq2.charAt(k) == '-';
							}
						}
						dist[i][j] = diff;
						nonzero = nonzero || diff != 0;
					}
				}
			}
			System.out.println(bl + " had nonzero: " + nonzero);
			
			//write distance
			BufferedWriter out = new BufferedWriter(new FileWriter(new File(
					DIR + bl + ".dist.txt")));
			//write header
			for(int i = 0; i < genes.size(); i++) {
				out.write("\t" + genes.get(i));
			}
			out.write("\n");
			for(int i = 0; i < genes.size(); i++) {
				out.write(genes.get(i));
				for(int j = 0; j < genes.size(); j++) {
					out.write("\t" + dist[i][j]);
				}
				out.write("\n");
			}
			out.close();
		}
	}
	

}
