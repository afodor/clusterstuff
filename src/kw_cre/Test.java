/*
 * generates the fasta files and scripts to align beta-lactamase genes
 * Klebsiella pneumonia only
 */
package kw_cre;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;

public class Test {
	public static String DIR = "/nobackup/afodor_research/kwinglee/cre/chs_v_cards/";
	public static String BLAST_HITS = "blastCardsProHomologGeneTables_pid80.0_len90.0_collapsed.txt";
	public static String FASTA_DIR = "/nobackup/afodor_research/kwinglee/cre/rbh/geneFastas/";

	public static void main(String[] args) throws Exception {
		String[] genes = new String[]{"KPC", "LEN", "OKP", "OXA", "SHV", "TEM"};//genes to analyze
		
		String outDir = DIR + "betaLactamaseAlignments/";
		
		//get map of gene name to line of hits
		BufferedReader hits = new BufferedReader(new FileReader(new File(
				DIR + BLAST_HITS)));
		HashMap<String, String> geneHits = new HashMap<String, String>();//map of gene to the corresponding line in the hits file
		for(String line = hits.readLine(); line != null; line = hits.readLine()) {
			String[] sp = line.split("\t");
			geneHits.put(sp[0], line);
		}
		hits.close();
		
		//generate fasta file for each gene
		for(String g : genes) {
			if(!geneHits.containsKey(g)) {
				System.err.println("Results file is missing gene " + g);
			} else {
				String line = geneHits.get(g);
				String[] sp = line.split("\t");
				HashSet<String> matches = new HashSet<String>(Arrays.asList(sp));
				HashSet<String> written = new HashSet<String>();
				matches.remove(g);//remove first column
				matches.remove("NA");//remove NAs
				for(String m : matches) {
					String[] msp = m.split(";");//some genomes have multiple genes; include both
					for(String s : msp) {
						if(s.contains("klebsiella_pneu")) {//kleb pneu only
							written.add(s);
						} else {
							System.out.println(s);
						}
					}
				}
				//get fasta genes
				BufferedReader fasta = new BufferedReader(new FileReader(new File(
						outDir + g + ".fasta")));
				HashSet<String> inFasta = new HashSet<String>();
				for(String line1 = fasta.readLine(); line1 != null; line1 = fasta.readLine()) {
					if(line1.startsWith(">")) {
						inFasta.add(line1.replace(">", ""));
					}
				}
				fasta.close();
				
				//compare sets
				System.out.println(g + "\thits: " + written.size() + "\tfasta: " + inFasta.size());
			}
		}
		
	}
}

