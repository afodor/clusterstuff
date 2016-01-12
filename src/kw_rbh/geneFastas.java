/**
 * For each Broad genome, generate a fasta file of all genes from the gtf file
 */
package kw_rbh;

import java.io.File;

public class geneFastas {
	public static String GenomeTopDir = "/nobackup/afodor_research/af_broad";

	public static void main(String[] args) {
		analyzeFolder("carolina");
		analyzeFolder("susceptible");
		analyzeFolder("resistant");
	}
	
	//analyze the genomes in the given folder
	public static void analyzeFolder(String folder) {
		String outDir = "/nobackup/afodor_research/kwinglee/cre/rbh/" + folder;
		File dir = new File(GenomeTopDir + folder);
		File[] allFiles = dir.listFiles();
		for(File gtf : allFiles) {
			if(gtf.getName().endsWith(".genes.gtf")) {
				String name = gtf.getName().replace(".genes.gtf", "");//genome to analyze
				//make directory for output files
				File gtfDir = new File(outDir + name);
				gtfDir.mkdirs();
				
				//get corresponding fasta file
				
				//make gene files
			}
		}
	}
}
