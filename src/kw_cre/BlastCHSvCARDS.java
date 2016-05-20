/*
 * blast genes and scaffolds of all CHS genomes against CARDS database
 */
package kw_cre;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;

public class BlastCHSvCARDS {
	public static String OUTDIR = "/nobackup/afodor_research/kwinglee/cre/chs_v_cards/";
	public static String SCAFFDIR = "/nobackup/afodor_research/af_broad/carolina/";
	public static String GENEDIR = "/nobackup/afodor_research/kwinglee/cre/rbh/carolina/";
	public static String CARDS = "/users/kwinglee/card/nucleotide_fasta.protein_homolog.fasta";
	
	public static void main(String[] args) throws IOException {
		//script for scaffolds
		String[] files = new File(SCAFFDIR).list();
		BufferedWriter script = new BufferedWriter(new FileWriter(new File(
				OUTDIR + "blastCHSscaffToCARDS")));
		script.write("module load blast\n");
		for(String f : files) {
			if(f.endsWith(".scaffolds.fasta")) {
				script.write("blastn -query " + SCAFFDIR + f 
						+ " -db " + CARDS +
						" -outfmt 7 -out " + OUTDIR + 
						f.replace(".fasta", "") + "_v_cardsProHomolog" + 
						"\n");
			}
		}		
		script.close();
		
		//script for genes
		files = new File(GENEDIR).list();
		script = new BufferedWriter(new FileWriter(new File(
				OUTDIR + "blastCHSgenesToCARDS")));
		script.write("module load blast\n");
		for(String f : files) {
			if(f.endsWith("_allGenes.fasta")) {
				script.write("blastn -query " + GENEDIR + f 
						+ " -db " + CARDS +
						" -outfmt 7 -out " + OUTDIR + 
						f.replace("_allGenes.fasta", ".genes") + "_v_cardsProHomolog" + 
						"\n");
			}
		}
		script.close();
	}

}
