/*
 * Generate scripts to run kraken 
 * 12/9/16
 */
package kw_machineLearning;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;

public class KrakenObesity {
	public static String KRAKEN_DIR = "/nobackup/afodor_research/kwinglee/software/kraken/";
	//public static String DB = KRAKEN_DIR + "krakenStandardDB";
	public static String DB = KRAKEN_DIR + "minikraken_20141208";
	public static String FASTQ_DIR = "/nobackup/afodor_research/kwinglee/machineLearning/obesity/fastqs/";
	public static String SCRIPT_DIR = "/projects/afodor_research/kwinglee/scripts/machineLearning/obesity/";
	//public static String OUT_DIR = "/nobackup/afodor_research/kwinglee/machineLearning/obesity/stdKraken/";
	public static String OUT_DIR = "/nobackup/afodor_research/kwinglee/machineLearning/obesity/minikraken/";
	
	public static void main(String[] args) throws IOException {
		File out = new File(OUT_DIR);
		if(!out.exists()) {
			out.mkdirs();
		}
		//set up script to run everything
		BufferedWriter script = new BufferedWriter(new FileWriter(new File(
				//SCRIPT_DIR + "stdKrakenObesity")));
				SCRIPT_DIR + "minikrakenObesity")));
		script.write("#PBS -l walltime=600:00:00\n");
		script.write("#PBS -l mem=10GB\n");
		String[] fastqs = new File(FASTQ_DIR).list();
		for(String fq : fastqs) {
			if(fq.endsWith(".1.fq.gz")) {
				String name = fq.replace(".1.fq.gz", "");
				//String seqName = OUT_DIR + name + "_stdKraken";
				String seqName = OUT_DIR + name + "_minikraken";
				
				//run kraken
				script.write(KRAKEN_DIR + "kraken --preload --db " 
						+ DB + " --fastq-input --gzip-compressed " +
						FASTQ_DIR + fq + " --threads 2 > " + seqName + "\n");
				//translate output
				script.write(KRAKEN_DIR + "kraken-translate --db "
						+ DB + " " + seqName + " > " + seqName + "_translate\n");
				script.write(KRAKEN_DIR + "kraken-translate --mpa-format --db "
						+ DB + " " + seqName + " > " + seqName + "_mpa\n");
			}
		}
		
		script.close();
	}

}
