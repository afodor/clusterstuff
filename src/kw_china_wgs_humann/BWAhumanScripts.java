/*
 * Scripts to use BWA to align China WGS reads against the human hg38 genome
 */
package kw_china_wgs_humann;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;

public class BWAhumanScripts {
	public static final String BASE_DIR = "/nobackup/afodor_research/kwinglee/china/wgs/";
	public static final String REF = "/nobackup/afodor_research/kwinglee/hg38/hg38.fa";
	
	public static void main(String[] args) throws IOException {
		String scriptDir = BASE_DIR + "bwaScripts/";//folder for scripts
		new File(scriptDir).mkdirs();
		String fastaDir = BASE_DIR + "fastas/";//folder with fasta files
		String outDir = BASE_DIR + "alignToHG38/";//folder to write blast results
		new File(outDir).mkdirs();
		
		//script to run all files
		BufferedWriter runAll = new BufferedWriter(new FileWriter(new File(
				scriptDir + "runAll.sh")));
		
		File[] fastas = new File(fastaDir).listFiles();
		for(File f : fastas) {
			if(f.getName().endsWith(".fa")) {
				String name = f.getName().replace(".fa", "");
				String scriptName = "align_" + name;
				
				//add to runAll
				runAll.write("qsub -q \"viper_batch\" " + scriptName + "\n");
				
				//write individual script
				BufferedWriter script = new BufferedWriter(new FileWriter(new File(
						scriptDir + scriptName)));
				script.write("#PBS -l walltime=500:00:00\n");
				script.write("module load bwa\n");
				script.write("bwa mem " + REF + f.getAbsolutePath());//command to align file
				script.close();
			}
		}
		
		runAll.close();
	}

}
