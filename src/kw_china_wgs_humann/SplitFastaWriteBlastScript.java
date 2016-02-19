/*
 * Take whole genome sequencing fastas, split, write scripts to run blast
 */
package kw_china_wgs_humann;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;

public class SplitFastaWriteBlastScript {
	private static final int NUM_READS = 100000;//max number of reads per file
	private static final String ROOT_DIR = "/nobackup/afodor_research/kwinglee/china/wgs/";
	private static final String FASTA_DIR = ROOT_DIR + "fastas/";
	private static final String SCRIPT_DIR = ROOT_DIR + "humScripts/";
	private static final String BLAST_DIR = ROOT_DIR + "kegg_split_blastx_results/";
	private static final String SPLIT_DIR = ROOT_DIR + "splitFastas/";
	private static final String DB = "/nobackup/afodor_research/kwinglee/china/wgs/allKeggPep/allKegg";
	
	public static void main(String[] args) throws IOException {
		String[] fastas = new File(FASTA_DIR).list();
		BufferedWriter runAll = new BufferedWriter(new FileWriter(new File(
				SCRIPT_DIR + "runAllSplit.sh")));
		for(String fa : fastas) {
			if(fa.endsWith(".fa")) {
				String name = fa.replace("_1.fa", "");
				
				int readCount = 0;
				int fileCount = 0;
				BufferedReader br = new BufferedReader(new FileReader(new File(
						FASTA_DIR + fa)));
				String head = br.readLine();
				//set up new fasta file
				String newFile = "split_" + name + "_" + fileCount;
				BufferedWriter fasta = new BufferedWriter(new FileWriter(new File(
						SPLIT_DIR + newFile + ".fa")));
				
				while(head != null) {
					if(readCount == NUM_READS) {//set up new file
						readCount = 0;
						//write blast script
						BufferedWriter script = new BufferedWriter(new FileWriter(new File(
								SCRIPT_DIR + "sBlast_" + newFile)));
						script.write("#PBS -l walltime=300:00:00\n");
						script.write("module load blast\n");
						script.write("blastx -outfmt 6 -db " + DB + " -query " +
								SPLIT_DIR + newFile + ".fa" + " -out " +
								BLAST_DIR + "kegg_" + newFile + ".txt\n");//blast command
						script.close();
						//add to run all
						runAll.write("qsub -q \"Cobra_batch\" " + "sBlast_" + newFile + "\n");

						//set up new fasta file
						fasta.close();
						fileCount++;
						newFile = "split_" + name.replace("_1", "") + "_" + fileCount;
						fasta = new BufferedWriter(new FileWriter(new File(
								SPLIT_DIR + newFile + "_" + fileCount + ".fa")));
					}
					//write new fasta files
					String seq = br.readLine();
					fasta.write(head + "\n" + seq + "\n");
					head = br.readLine();
					readCount++;
				}
				fasta.close();
				br.close();
				
				//write blast script
				BufferedWriter script = new BufferedWriter(new FileWriter(new File(
						SCRIPT_DIR + "sBlast_" + newFile)));
				script.write("#PBS -l walltime=300:00:00\n");
				script.write("module load blast\n");
				script.write("blastx -outfmt 6 -db " + DB + " -query " +
						SPLIT_DIR + newFile + ".fa" + " -out " +
						BLAST_DIR + "kegg_" + newFile + ".txt\n");//blast command
				script.close();
				//add to run all
				runAll.write("qsub -q \"Cobra_batch\" sBlast_" + newFile + "\n");
				
			}
		}
		runAll.close();
	}

}
