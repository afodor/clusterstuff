/*
 * generates the scripts for stitching reads using fastq-join
 */
package kw_jobinDolphin;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;

public class QiimeReadStitchingScripts {
	public static String BASE_DIR = "/nobackup/afodor_research/kwinglee/jobin/dolphin/";
	public static String FQ_DIR = BASE_DIR + "trimmed_fastqs/";
	public static String OUT_DIR = BASE_DIR + "stitched_reads/";
	public static String SCRIPT_DIR = BASE_DIR + "qiimeScripts/";
	
	public static void main(String[] args) throws IOException {
		//get list of files
		String[] fqs = new File(FQ_DIR).list();
		
		//set up runAll
		BufferedWriter runAll = new BufferedWriter(new FileWriter(new File(
				SCRIPT_DIR + "allStitch.sh")));
		
		for(String fq: fqs) {
			if(fq.endsWith("_R1.fastq")) {
				String sample = fq.replace("_R1.fastq", "");
				String name = "stitch_" + sample;
				//individual script
				BufferedWriter script = new BufferedWriter(new FileWriter(new File(
						SCRIPT_DIR + name)));
				/*//for qiime
				script.write("module load openmpi\n");
				script.write("module load qiime\n");
				script.write("join_paired_ends.py -f " + FQ_DIR + fq
						+ " -r " + FQ_DIR + sample + "_R2.fastq -o " +
						OUT_DIR + sample + "%.fastq\n");*/
				script.write("/users/kwinglee/ea-utils-read-only/clipper/fastq-join " +
						FQ_DIR + fq + " " + FQ_DIR + sample + "_R2.fastq" + 
						" -o " + OUT_DIR + sample  + "." + "%.fastq");
				script.close();
				
				//add to runAll
				runAll.write("qsub -q \"Cobra\" " + name + "\n");
			}
		}
		runAll.close();
	}

}
