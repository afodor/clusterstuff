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

public class AlignSelectCardsHits {
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
				matches.remove(g);//remove first column
				matches.remove("NA");//remove NAs
				BufferedWriter fasta = new BufferedWriter(new FileWriter(new File(
						outDir + g + ".fasta")));
				for(String m : matches) {
					String[] msp = m.split(";");//some genomes have multiple genes; include both
					for(String s : msp) {
						if(s.contains("klebsiella_pneu")) {//kleb pneu only
							//get fasta sequence
							String[] name = s.split("_A[E|F]");
							if(name.length != 2) {
								fasta.close();
								throw new Exception("bad gene split " + name.length + " " + s);
							}
							BufferedReader seq = new BufferedReader(new FileReader(new File(
									FASTA_DIR + name[0] + File.separator + s + ".fasta")));
							for(String seqline = seq.readLine(); seqline != null; seqline = seq.readLine()) {
								fasta.write(seqline + "\n");
							}
							seq.close();
						}
					}
				}
				fasta.close();
			}
		}
		
		//generate scripts
		BufferedWriter script = new BufferedWriter(new FileWriter(new File(
				outDir + "alignCommands")));
		script.write("module load mafft\n");
		for(String g : genes) {
			script.write("echo \"\"\n");
			script.write("echo \"Starting " + g + "\"\n");
			script.write("echo \"\"\n");
			script.write("mafft --auto --adjustdirection " + outDir + g + ".fasta" 
						+ " > " + outDir + g + ".mafftAlign.fast\n");
		}
		script.close();
	}
}

