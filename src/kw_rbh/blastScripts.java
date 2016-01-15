/**
 * code to generate scripts to blast each genome against each other genome
 */
package kw_rbh;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;

public class blastScripts {
	public static String DIR = "/nobackup/afodor_research/kwinglee/cre/rbh/";
	
	public static void main(String[] args) throws IOException {
		BufferedWriter runAll = new BufferedWriter(new FileWriter(new File(
				DIR + "makedbScripts/runAll.sh")));//script to run all files
		String[] folders = {"carolina", "susceptible", "resistant"};
		for(String f1 : folders) {
			for(String f2 : folders) {
				//folder to put comparisons in
				String resultsFolder = DIR + "blastResults/" + f1 + "_v_" + f2;
				File resFile = new File(resultsFolder);
				resFile.mkdirs();
				
				//make script that runs each set of comparisons separately
				BufferedWriter compare = new BufferedWriter(new FileWriter(new File(
						DIR + "blastScripts/run_" + f1 + "_v_" + f2 + ".sh")));
				runAll.write("./run_" + f1 + "_v_" + f2 + ".sh");
				
				//scripts to compare each genome against each other genome
				File[] genomes1 = new File(DIR + f1).listFiles();
				File[] genomes2 = new File(DIR + f2).listFiles();
				for(File g1 : genomes1) {
					for(File g2 : genomes2) {
						String gen1 = g1.getName();
						String gen2 = g2.getName();
						
						//set up individual script
						BufferedWriter script = new BufferedWriter(new FileWriter(new File(
								DIR + "makedbScripts/blast_" + gen1 + "_v_" + gen2)));
						script.write("module load blast\n");
						script.write("blastn -query " + g1.getAbsolutePath() + "/" + gen1 + "_allGenes.fasta -db " + 
								g2.getAbsolutePath() + "/" + gen2 + "_allGenes.fasta -outfmt 7 -out " +
								resultsFolder + "/" + gen1 + "_v_" + gen2 + ".txt");
						script.close();
						
						//add to runAll
						runAll.write("qsub -q \"viper_batch\" blast_" + gen1 + "_v_" + gen2 + "\n");
					}
				}
				compare.close();
			}
		}
		runAll.close();
	}

}
