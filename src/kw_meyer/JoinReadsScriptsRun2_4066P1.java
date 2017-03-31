/*
 * generate the scripts to join/stitch the run2 4066 p1 reads
 */
package kw_meyer;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.HashSet;

public class JoinReadsScriptsRun2_4066P1 {
	//directory containing the needed files; also where will write results
	public static final String DIR = "/nobackup/afodor_research/kwinglee/meyer/";
	public static final String SCRIPTDIR = "/projects/afodor_research/kwinglee/scripts/meyer/run2/joinScripts/";
	
	public static void main(String[] args) throws IOException {
		//get list of samples to join
		HashSet<String> newSeqs = new HashSet<String>();
		File[] samps = new File(DIR + "cardia_seq2_fastqs/4066-KAM-P1-34087111/").listFiles();
		for(File f : samps) {
			String[] reads = f.list();
			for(String r : reads) {
				newSeqs.add(r.replace(".gz", ""));
			}
		}
		
		//write scripts
		File[] fastqs = new File(DIR + "run2filteredSeqs").listFiles();
		BufferedWriter allScript = new BufferedWriter(new FileWriter(new File(
				SCRIPTDIR + "joinAll4066P1.sh")));
		
		for(File f : fastqs) {
			if(f.getName().endsWith("_R1_001.fastq") && newSeqs.contains(f.getName())) {
				String sampleID = f.getName().replace("_L001_R1_001.fastq", "");
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
