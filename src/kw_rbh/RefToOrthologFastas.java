/*
 * Make fasta files of all orthologs for each gene in CHS11 if that gene
 * has more orthologs than the cutoff 
 */
package kw_rbh;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

public class RefToOrthologFastas {
	public static final int CUTOFF = 200;//(determined from ortholog_cutoff_plot.R)
	public static String DIR = "/nobackup/afodor_research/kwinglee/cre/rbh/rbhOrthologs/";

	public static void main(String[] args) throws IOException {
		////get set of orthologs for each gene
		BufferedReader br = new BufferedReader(new FileReader(new File(
				"orthologTables/orthologNameTable_carolina_klebsiella_pneumoniae_chs_11.0.txt")));
		String line = br.readLine();

		//initialize the sets with the genes
		List<Set<String>> sets = new ArrayList<Set<String>>();
		String[] genes = line.split("\t");
		for(int i = 1; i < genes.length; i++) {//skip first element, aka the genome ids
			Set<String> s = new HashSet<String>();
			s.add(genes[i]);
			sets.add(s);
		}
		line = br.readLine();

		//add genes to set
		while(line != null) {
			String[] sp = line.split("\t");
			for(int i = 1; i < sp.length; i++) {
				sets.get(i-1).add(sp[i]);
			}
			line = br.readLine();
		}
		br.close();

		//remove NA
		for(int i = 0; i < sets.size(); i++) {
			sets.get(i).remove("NA");
		}

		////if set size is greater than cutoff, print all genes to a fasta file
		String fastaDir = DIR + "chs11OrthogroupFastas/";
		new File(fastaDir).mkdirs();
		for(int i=0; i < sets.size(); i++) {
			Set<String> set = sets.get(i);
			if(set.size() > CUTOFF) {
				BufferedWriter fasta = new BufferedWriter(new FileWriter(new File(
						fastaDir + "orthogroups_" + genes[i+1] + ".fasta")));
				for(String s : set) {
					//String genome = s.replaceFirst("_[A-Z]*[0-9]*_[0-9]*$", "");//works on carolina, haven't checked other groups
					String[] sp = s.split("_");
					String genome = "";
					for(int j = 0; j < sp.length-2; j++) {
						genome += sp[j] + "_";
					}
					genome = genome.replaceFirst("_$", "");
					br = new BufferedReader(new FileReader(new File(
							"/nobackup/afodor_research/kwinglee/cre/rbh/geneFastas/" + 
									genome + "/" + s + ".fasta")));
					line = br.readLine();
					while(line != null) {
						fasta.write(line + "\n");
						line = br.readLine();
					}
					br.close();
				}
				fasta.close();
			}
		}
	}
}
