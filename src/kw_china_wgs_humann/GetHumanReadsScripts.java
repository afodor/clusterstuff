/*
 * Scripts to get list of China WGS reads that mapped to human genome
 * Also generates stats tables
 */
package kw_china_wgs_humann;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;

public class GetHumanReadsScripts {
	public static final String BASE_DIR = "/nobackup/afodor_research/kwinglee/china/wgs/";
	
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
				String scriptName = "filter_" + name;
				
				//add to runAll
				runAll.write("qsub -q \"viper_batch\" " + scriptName + "\n");
				
				//write individual script
				BufferedWriter script = new BufferedWriter(new FileWriter(new File(
						scriptDir + scriptName)));
				script.write("module load samtools\n");
				script.write("samtools view -F 4 " + outDir + name + ".sam > " 
						+ outDir + name + ".human.mapped.sam"); //get mapped reads; use -f 4 for unmapped
				script.write("samtools flagstat " + outDir + name + ".sam > "
						+ outDir + name + ".stats.txt");//get stats
				script.close();
			}
		}
		
		runAll.close();
	}

}
