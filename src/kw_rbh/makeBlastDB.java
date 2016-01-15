/**
 * generates the scripts to make nucleotide blast databases from the files containing all genes
 */

package kw_rbh;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;

public class makeBlastDB {
	public static String DIR = "/nobackup/afodor_research/kwinglee/cre/rbh/";

	public static void main(String[] args) throws IOException {
		BufferedWriter runAll = new BufferedWriter(new FileWriter(new File(
				DIR + "makedbScripts/runAll.sh")));//script to run all files
		String[] folders = {"carolina", "susceptible", "resistant"};
		/*for(String f : folders) {
			File folder = new File(DIR + f);
			File[] genomes = folder.listFiles();
			for(File g : genomes) {
				String gen = g.getName();
				
				//set up individual script
				BufferedWriter script = new BufferedWriter(new FileWriter(new File(
						DIR + "makedbScripts/run" + gen)));
				script.write("module load blast\n");
				script.write("makeblastdb -in " + g.getAbsolutePath() + "/" + f + "_" + gen + "_allGenes.fasta -dbtype 'nucl'\n");
				script.close();
				
				//add to runAll
				runAll.write("qsub -q \"viper_batch\" run" + gen + "\n");
			}
		}*/
		for(String f : folders) {
			File folder = new File(DIR + f);
			File[] genomes = folder.listFiles();
			for(File g : genomes) {
				String gen = g.getName();
				
				if(gen.endsWith("allGenes.fasta")) {
					//set up individual script
					BufferedWriter script = new BufferedWriter(new FileWriter(new File(
							DIR + "makedbScripts/run" + gen)));
					script.write("module load blast\n");
					/*script.write("makeblastdb -in " + g.getAbsolutePath() + "/" + f + "_" + gen + "_allGenes.fasta -dbtype 'nucl'\n");*/
					script.write("makeblastdb -in " + g.getAbsolutePath() + "/" + gen + " -dbtype 'nucl'\n");
					script.close();
					
					//add to runAll
					runAll.write("qsub -q \"viper_batch\" run" + gen + "\n");
				}
			}
		}
		runAll.close();
	}
}
