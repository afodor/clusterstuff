/*
 * generate the scripts to join/stitch the demultiplxed reads
 */
package kw_jobinApcTumor;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;

public class JoinReadsScripts {
	//directory containing the needed files; also where will write results
	public static final String DIR = "/nobackup/afodor_research/kwinglee/jobin/apcTumor/";
	public static final String SCRIPTDIR = "/projects/afodor_research/kwinglee/scripts/jobin/apcTumor/joinScripts/";
	
	public static void main(String[] args) throws IOException {
		File[] fastqs = new File(DIR + "fastqs").listFiles();
		BufferedWriter allScript = new BufferedWriter(new FileWriter(new File(
				SCRIPTDIR + "joinAll.sh")));
		for(File f : fastqs) {
			if(f.getName().endsWith("_R1.fastq") && !f.getName().startsWith("other")) {
				String sampleID = f.getName().replace("_R1.fastq", "");
				String scriptName = "join_" + sampleID;
				
				String r2 = f.getAbsolutePath().replace("_R1.fastq", "_R2.fastq");
				BufferedWriter script = new BufferedWriter(new FileWriter(new File(
						SCRIPTDIR + scriptName)));
				script.write("/users/kwinglee/ea-utils-read-only/clipper/fastq-join " +
						f.getAbsolutePath() + " " + r2 + 
						" -o " + DIR + "joinedReads/" + sampleID + "%.fastq");
				script.close();
				
				allScript.write("qsub -q \"copperhead\" " + scriptName + "\n");
			}
		}
		allScript.close();
	}
}
