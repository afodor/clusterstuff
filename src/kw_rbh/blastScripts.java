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
	public static int SUB_SIZE = 1000;//number of commands per subset file
	
	public static void main(String[] args) throws IOException {
		new File(DIR + "blastScripts/subsets").mkdirs();//folder to put files to run subsets
		BufferedWriter runAll = new BufferedWriter(new FileWriter(new File(
				DIR + "blastScripts/runAll.sh")));//script to run all files
		String[] folders = {"carolina", "susceptible", "resistant"};
		int numCmd = 0; //number of commands so far
		BufferedWriter subset = null; //writer to write subset of commands
		String subset_name = DIR + "blastScripts/subsets/subset_";
		for(String f1 : folders) {
			for(String f2 : folders) {
				//folder to put comparisons in
				String resultsFolder = DIR + "blastResults/" + f1 + "_v_" + f2;
				File resFile = new File(resultsFolder);
				resFile.mkdirs();
				
				//make script that runs each set of comparisons separately
				BufferedWriter compare = new BufferedWriter(new FileWriter(new File(
						DIR + "blastScripts/run_" + f1 + "_v_" + f2 + ".sh")));
				runAll.write("./run_" + f1 + "_v_" + f2 + ".sh\n");
				
				//scripts to compare each genome against each other genome
				File[] genomes1 = new File(DIR + f1).listFiles();
				File[] genomes2 = new File(DIR + f2).listFiles();
				for(File g1 : genomes1) {
					if(g1.getName().endsWith(".fasta")) {
						String gen1 = g1.getName().replace("_allGenes.fasta", "");
						//put all comparisons to this genome in their own folder (to make file system easier to access)
						File genResults = new File(resultsFolder + "/" + gen1);
						genResults.mkdirs();
						for(File g2 : genomes2) {
							if(g2.getName().endsWith(".fasta")) {
								if(numCmd % SUB_SIZE == 0) {
									if(numCmd != 0) {
										subset.close();
									} 
									subset = new BufferedWriter(new FileWriter(new File(
											subset_name + numCmd/SUB_SIZE + ".sh")));
								}
								
								String gen2 = g2.getName().replace("_allGenes.fasta", "");
								
								//set up individual script
								BufferedWriter script = new BufferedWriter(new FileWriter(new File(
										DIR + "blastScripts/blast_" + gen1 + "_v_" + gen2)));
								script.write("module load blast\n");
								script.write("blastn -query " + g1.getAbsolutePath() + " -db " + 
										g2.getAbsolutePath() + " -outfmt 7 -out " +
										genResults.getAbsolutePath() + "/" + gen1 + "_v_" + gen2 + ".txt\n");
								script.close();
								
								//add to runAll
								compare.write("qsub -q \"viper_batch\" blast_" + gen1 + "_v_" + gen2 + "\n");
								subset.write("qsub -q \"viper_batch\" blast_" + gen1 + "_v_" + gen2 + "\n");
							}
						}
					}
				}
				compare.close();
			}
		}
		subset.close();
		runAll.close();
		System.out.println("Total commands = " + numCmd);
	}

}
