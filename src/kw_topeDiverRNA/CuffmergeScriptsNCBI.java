/*
 * make cuffmerge input tables and scripts
 */
package kw_topeDiverRNA;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;

public class CuffmergeScriptsNCBI {
	public static String SCRIPTDIR = "/projects/afodor_research/kwinglee/scripts/tope/diverRNA/cufflinks/";
	private static String CUFFDIR = "/nobackup/afodor_research/kwinglee/tope/diverRNAseq/cufflinksHumanResultsNCBI/";
	private static String HUMGFF = "/nobackup/afodor_research/kwinglee/hg38/GCF_000001405.36_GRCh38.p10_genomic.gff";
	private static String HUMFA = "/nobackup/afodor_research/kwinglee/hg38/GCF_000001405.36_GRCh38.p10_genomic.fna";

	public static void main(String[] args) throws IOException {
		//set up manifest file
		String gtfList = CUFFDIR + "cuffmerge_assembly_list.txt";
		BufferedWriter out = new BufferedWriter(new FileWriter(new File(gtfList)));
		//make manifest files
		File[] files = new File(CUFFDIR).listFiles();
		for(File f : files) {
			if(f.isDirectory()) {
				out.write(f.getAbsolutePath() + File.separator + "transcripts.gtf\n");
			}
		}
		out.close();

		//set up scripts (for Cobra)
		BufferedWriter script = new BufferedWriter(new FileWriter(new File(
				SCRIPTDIR + "cuffmergeNCBI")));
		script.write("#PBS -l procs=1\n");
		script.write("PATH=$PATH:/nobackup/afodor_research/kwinglee/software/cufflinks-2.2.1.Linux_x86_64\n");
		script.write("cuffmerge -o " + CUFFDIR + "cuffmergeNCBI" + " -g " + HUMGFF
				+ " -p 4 -s " + HUMFA + " " + gtfList + "\n");
		script.close();
	}
}
