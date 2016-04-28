/*
 * Generate scripts to run kraken using the minikraken database
 * 4/28/17
 */
package kw_china_wgs_kraken;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;

public class MinikrakenScripts {
	public static String KRAKEN_DIR = "/nobackup/afodor_research/kwinglee/software/kraken/";
	public static String DB = KRAKEN_DIR + "minikraken_20141208";
	public static String FASTA_DIR = "/nobackup/afodor_research/kwinglee/china/wgs/fastas/";
	public static String SCRIPT_DIR = "/nobackup/afodor_research/kwinglee/china/wgs/krakenScripts/";
	public static String OUT_DIR = "/nobackup/afodor_research/kwinglee/china/wgs/minikrakenResults";
	
	public static void main(String[] args) throws IOException {
		//set up script to run everything
		BufferedWriter runAll = new BufferedWriter(new FileWriter(new File(
				SCRIPT_DIR + "runMiniKraken.sh")));
		String[] fastas = new File(FASTA_DIR).list();
		for(String fa : fastas) {
			if(fa.endsWith(".fa")) {
				String name = fa.replace("_1.fa", "");
				String scriptName = "mini_" + name;
				String seqName = OUT_DIR + "minikrakenSeqs_" + name;
				
				//generate script
				BufferedWriter script = new BufferedWriter(new FileWriter(new File(
						SCRIPT_DIR + scriptName)));
				//run kraken
				script.write(KRAKEN_DIR + "kraken --preload --db " 
						+ DB + " --fasta-input " +
						FASTA_DIR + fa + " > " + seqName + "\n");
				//translate output
				script.write(KRAKEN_DIR + "kraken-translate --db "
						+ DB + " " + seqName + " > " + seqName + "_translate\n");
				script.write(KRAKEN_DIR + "kraken-translate --mpa-format --db "
						+ DB + " " + seqName + " > " + seqName + "_mpa\n");
				script.close();
				
				//add to runAll
				runAll.write("qsub -q \"Cobra\" " + scriptName + "\n");
			}
		}
		
		runAll.close();
	}

}
