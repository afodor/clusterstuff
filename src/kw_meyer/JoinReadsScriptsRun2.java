/*
 * generate the scripts to join/stitch the run2 reads
 */
package kw_meyer;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;

public class JoinReadsScriptsRun2 {
	//directory containing the needed files; also where will write results
	public static final String DIR = "/nobackup/afodor_research/kwinglee/meyer/";
	public static final String SCRIPTDIR = "/projects/afodor_research/kwinglee/scripts/meyer/run2/joinScripts/";
	
	public static void main(String[] args) throws IOException {
		File[] fastqs = new File(DIR + "run2filteredSeqs").listFiles();
		BufferedWriter allScript = new BufferedWriter(new FileWriter(new File(
				SCRIPTDIR + "joinAll.sh")));
		for(File f : fastqs) {
			if(f.getName().endsWith("_R1_001.fastq")) {
				String sampleID = f.getName().split("_")[0];
				String scriptName = "join_" + sampleID;
				
				String r2 = f.getAbsolutePath().replace("_R1_001.fastq", "_R2_001.fastq");
				BufferedWriter script = new BufferedWriter(new FileWriter(new File(
						SCRIPTDIR + scriptName)));
				script.write("#PBS -l procs=1\n");
				script.write("/users/kwinglee/ea-utils-read-only/clipper/fastq-join " +
						f.getAbsolutePath() + " " + r2 + 
						" -o " + DIR + "run2joinedReads/" + sampleID + "%.fastq\n");
				script.close();
				
				allScript.write("qsub -q \"copperhead\" " + scriptName + "\n");
			}
		}
		allScript.close();
	}
}
