/*
 * Take whole mouse and rRNA filtered fastas, split, write scripts to run blast
 */
package kw_jobinBiofilm_rnaseq;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;

public class SplitFastaWriteBlastScript {
	private static final int NUM_READS = 100000;//max number of reads per file
	private static final int NUM_JOBS = 800;//max number of jobs per runAll
	private static final String ROOT_DIR = "/nobackup/afodor_research/kwinglee/jobin/biofilm/rnaseq/";
	private static final String FASTA_DIR = ROOT_DIR + "mouseAndSilvaFiltered/";
	private static final String SCRIPT_DIR = "/projects/afodor_research/kwinglee/scripts/jobin/biofilmRNAseq/humScripts/";
	private static final String BLAST_DIR = ROOT_DIR + "kegg_split_blastx_mouseSilvaFiltered/";
	private static final String SPLIT_DIR = ROOT_DIR + "splitMouseSilvaFilteredFastas/";
	private static final String DB = "/nobackup/afodor_research/kwinglee/china/wgs/allKeggPep/allKegg";
	
	public static void main(String[] args) throws IOException {
		String[] fastas = new File(FASTA_DIR).list();
		int numJobs = 0;
		int numAll = 0;
		BufferedWriter runAll = new BufferedWriter(new FileWriter(new File(
				SCRIPT_DIR + "runAllSplit" + numAll + ".sh")));
		for(String fa : fastas) {
			if(fa.endsWith(".fasta")) {
				String fastaName = fa.replace(".mouseFiltered.silvaFiltered.fasta", "");
				
				int readCount = 0;
				int fileCount = 0;
				BufferedReader br = new BufferedReader(new FileReader(new File(
						FASTA_DIR + fa)));
				String head = br.readLine();
				//set up new fasta file
				String newFile = fastaName + "_split" + fileCount;
				BufferedWriter fasta = new BufferedWriter(new FileWriter(new File(
						SPLIT_DIR + newFile + ".fa")));
				
				while(head != null) {
					if(readCount == NUM_READS) {//set up new file
						readCount = 0;
						//write blast script
						String scriptName = "blastx_" + newFile;
						BufferedWriter script = new BufferedWriter(new FileWriter(new File(
								SCRIPT_DIR + scriptName)));
						script.write("#PBS -l walltime=400:00:00\n");
						//script.write("#PBS -l mem=30GB\n");
						script.write("module load blast\n");
						script.write("blastx -outfmt 6 -db " + DB + " -query " +
								SPLIT_DIR + newFile + ".fa" + " -out " +
								BLAST_DIR + "kegg_" + newFile + ".txt\n");//blast command
						script.close();
						//add to run all
						if(numJobs == NUM_JOBS) {
							runAll.close();
							numAll++;
							runAll = new BufferedWriter(new FileWriter(new File(
									SCRIPT_DIR + "runAllSplit" + numAll + ".sh")));
							numJobs = 0;
						}
						numJobs++;
						runAll.write("qsub -q \"Cobra_batch\" " + scriptName + "\n");

						//set up new fasta file
						fasta.close();
						fileCount++;
						newFile = fastaName + "_split" + fileCount;
						fasta = new BufferedWriter(new FileWriter(new File(
								SPLIT_DIR + newFile + ".fa")));
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
				String scriptName = "blastx_" + newFile;
				BufferedWriter script = new BufferedWriter(new FileWriter(new File(
						SCRIPT_DIR + scriptName)));
				script.write("#PBS -l walltime=400:00:00\n");
				//script.write("#PBS -l mem=30GB\n");
				script.write("module load blast\n");
				script.write("blastx -outfmt 6 -db " + DB + " -query " +
						SPLIT_DIR + newFile + ".fa" + " -out " +
						BLAST_DIR + "kegg_" + newFile + ".txt\n");//blast command
				script.close();
				//add to run all
				runAll.write("qsub -q \"Cobra_batch\" " + scriptName + "\n");
				
			}
		}
		runAll.close();
		System.out.println("Number scripts: " + numAll);
	}

}
