/**
 * For the given genome, print a table of scaffold, start and stop for each gene
 * Also get sizes of each scaffold
 * For graphs of ortholog p-values
 */
package kw_rbh;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;

public class GenomePositions {
	public static String GenomeDir = "/nobackup/afodor_research/af_broad/";
	public static String OutDir = "/nobackup/afodor_research/kwinglee/cre/rbh/";
	public static String Group = "carolina";
	public static String Genome = "klebsiella_pneumoniae_chs_11.0";//CHS11
	
	public static void main(String[] args) throws IOException {
		//get table of gene positions
		BufferedReader gtf = new BufferedReader(new FileReader(new File(
				GenomeDir + Group + "/" + Genome + ".genes.gtf")));
		BufferedWriter out = new BufferedWriter(new FileWriter(new File(
				OutDir + Group + "_" + Genome + "_genePositions.txt")));
		out.write("geneID\tscaffold\tgeneStart\tgeneStop\n");
		String line = gtf.readLine();
		while(line != null) {
			String[] sp = line.split("\t");
			if(sp[2].equals("exon")) { //only include exons (which include start and stop codon)
				//use gene_id
				String id = sp[8].split(";")[0].split("\"")[1];
				out.write(Group + "_" + Genome + "_" + id + "\t"
						+ sp[0] + "\t" + sp[3] + "\t" + sp[4] + "\n");
			}
			line = gtf.readLine();
		}
		
		gtf.close();
		out.close();
		
		//get scaffold sizes
		BufferedReader fasta = new BufferedReader(new FileReader(new File(
				GenomeDir + Group + "/" + Genome + ".scaffolds.fasta")));
		out = new BufferedWriter(new FileWriter(new File(
				OutDir + Group + "_" + Genome + "_scaffoldSizes.txt")));
		out.write("scaffold\tlength\n");
		line = fasta.readLine();
		int len = 0;//length of the sequence
		String scaff = "";
		while(line != null) {
			if(line.startsWith(">")) {
				if(len != 0) {
					out.write(scaff + "\t" + len + "\n");
				}
				scaff = line.replace(">", "");
				len = 0;
			} else {
				len += line.length();
			}
			line = fasta.readLine();
		}
		out.write(scaff + "\t" + len + "\n");
		fasta.close();
		out.close();
	}

}
