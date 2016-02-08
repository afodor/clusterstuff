/*
 * Scripts to run blastx on the China WGS reads against the KEGG database
 */
package kw_china_wgs_humann;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;

public class KEGGblastxScripts {
	public static final String BASE_DIR = "/nobackup/afodor_research/kwinglee/china/wgs/";
	//public static final String DB = "/nobackup/afodor_research/rbarner/kegg/allKegg"; 
	//public String DB = "/nobackup/afodor_research/kegg/seq_pep/seq_pep_database/kegg_pep.fasta";
	public static final String DB = "/nobackup/afodor_research/kwinglee/china/wgs/allKeggPep/allKegg";
	
	public static void main(String[] args) throws IOException {
		String scriptDir = BASE_DIR + "humScripts/";//folder for scripts
		new File(scriptDir).mkdirs();
		String fastaDir = BASE_DIR + "fastas/";//folder with fasta files
		String outDir = BASE_DIR + "kegg_blastx_results/";//folder to write blast results
		new File(outDir).mkdirs();
		
		//script to run all files
		BufferedWriter runAll = new BufferedWriter(new FileWriter(new File(
				scriptDir + "runAll.sh")));
		
		File[] fastas = new File(fastaDir).listFiles();
		for(File f : fastas) {
			if(f.getName().endsWith(".fa.gz")) {
				String name = f.getName().replace(".fa.gz", "");
				String scriptName = "blast_" + name;
				
				//add to runAll
				runAll.write(scriptName + "\n");
				
				//write individual script
				BufferedWriter script = new BufferedWriter(new FileWriter(new File(
						scriptDir + scriptName)));
				script.write("#PBS -l walltime=300:00:00\n");
				script.write("module load blast\n");
				script.write("gunzip -c " + f.getAbsolutePath() + " > " + 
						f.getAbsolutePath().replace(".gz", "") + "\n");//unzip fasta file
				script.write("blastx -outfmt 6 -db " + DB + " -query " +
						f.getAbsolutePath() + " -out " +
						outDir + "kegg_" + name + ".txt\n");//blast command
				script.close();
			}
		}
		
		runAll.close();
	}

}
