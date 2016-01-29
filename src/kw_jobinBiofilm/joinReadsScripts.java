/*
 * generate the scripts to join/stitch the demultiplxed reads
 */
package kw_jobinBiofilm;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;

public class joinReadsScripts {
	//directory containing the needed files; also where will write results
	public static final String DIR = "/nobackup/afodor_research/kwinglee/jobin/biofilm/";
	
	public static void main(String[] args) throws IOException {
		File[] fastqs = new File(DIR + "fastqs").listFiles();
		String scriptDir = DIR + "qiimeScripts/";
		BufferedWriter allScript = new BufferedWriter(new FileWriter(new File(
				scriptDir + "joinAll.sh")));
		for(File f : fastqs) {
			if(f.getName().endsWith("_R1.fastq")) {
				String sampleID = f.getName().split("_")[0];
				String scriptName = "join_" + sampleID;
				
				String r2 = f.getAbsolutePath().replace("_R1.fastq", "_R2.fastq");
				BufferedWriter script = new BufferedWriter(new FileWriter(new File(
						scriptDir + scriptName)));
				script.write("/users/kwinglee/ea-utils-read-only/clipper/fastq-join " +
						f.getAbsolutePath() + " " + r2 + 
						" -o " + DIR + "join/biofilm_" + sampleID + "%.fastq");
				script.close();
				
				allScript.write("qsub -q \"viper_batch\" " + scriptName + "\n");
			}
		}
		allScript.close();
	}
}
