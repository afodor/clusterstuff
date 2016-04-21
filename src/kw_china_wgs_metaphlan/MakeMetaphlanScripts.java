/*
 * create scripts to run metaphlan for all China WGS
 */
package kw_china_wgs_metaphlan;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;

public class MakeMetaphlanScripts {
	public static String BASE_DIR = "/nobackup/afodor_research/kwinglee/china/wgs/";//directory for outputs
	public static String MET_DIR = "/nobackup/afodor_research/kwinglee/software/metaphlan2";//location of metaphlan
	public static String BOW_DIR = "/nobackup/afodor_research/kwinglee/software/bowtie2-2.2.8";//bowtie2 location
	
	public static void main(String[] args) throws IOException {
		String scriptDir = BASE_DIR + "metScripts/";
		String outDir = BASE_DIR + "metaphlanResults/";
		String fastaDir = BASE_DIR + "fastas/";
		
		BufferedWriter runAll = new BufferedWriter(new FileWriter(new File(
				scriptDir + "runAll.sh")));
		String[] fastas = new File(fastaDir).list();
		for(String f : fastas) {
			if(f.endsWith(".fa")) {
				String sample = f.replace("_1.fa", "");
				String scriptName = "met_" + sample;
				BufferedWriter script = new BufferedWriter(new FileWriter(new File(
						scriptDir + scriptName)));
				//load python
				script.write("module load python/2.7.10\n");
				//set path variables
				script.write("export PATH=" + MET_DIR + ":" + BOW_DIR + ":$PATH\n");
				script.write("export mpa_dir=" + MET_DIR + "\n");
				//run metaphlan
				script.write("metaphlan2.py " + fastaDir + f + 
						" --bowtie2out " + outDir + "metaphlan_bowtie2_" + sample + ".bz2" +
						" --input_type fasta --nproc 2 > " + 
						outDir + "metaphlan_table_" + sample + "\n");
				script.close();
				
				runAll.write("qsub -q \"Cobra\" " + scriptName + "\n");
			}
		}
		
		runAll.close();
	}

}
