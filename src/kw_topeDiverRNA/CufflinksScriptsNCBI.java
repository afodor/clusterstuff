/*
 * generate the scripts to run cufflinks on Tophat results
 */
package kw_topeDiverRNA;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;

public class CufflinksScriptsNCBI {
	private static String BASEDIR = "/nobackup/afodor_research/kwinglee/tope/diverRNAseq/";
	private static String TOPHATDIR = BASEDIR + "tophatAlignToHumanNCBI/";
	private static String CUFF = "/nobackup/afodor_research/kwinglee/software/cufflinks-2.2.1.Linux_x86_64/cufflinks";
	private static String OUTDIR = BASEDIR + "cufflinksHumanResultsNCBI/";
	private static String SCRIPTDIR = "/projects/afodor_research/kwinglee/scripts/tope/diverRNA/cufflinks/";
	private static String GFF = "/nobackup/afodor_research/kwinglee/hg38/GCF_000001405.36_GRCh38.p10_genomic.gff";

	public static void main(String[] args) throws IOException {
		String[] sams = new File(TOPHATDIR).list();
		BufferedWriter runAll = new BufferedWriter(new FileWriter(new File(
				SCRIPTDIR + "runCufflinksTophatNCBI.sh")));
		for(String s : sams) {
			String scriptName = "ncbiCufflinks_" + s; 
			BufferedWriter script = new BufferedWriter(new FileWriter(new File(
					SCRIPTDIR + scriptName)));
			script.write("#PBS -l procs=1\n");
			script.write("#PBS -l walltime=500:00:00,mem=20GB\n");
			script.write(CUFF + " --GTF-guide " + GFF + " --output-dir " + 
					OUTDIR + s + " " + TOPHATDIR + s + 
					File.separator + "accepted_hits.bam\n");
			script.close();

			runAll.write("qsub -q \"copperhead\" " + scriptName + "\n");
			//runAll.write("qsub -q \"Cobra\" " + scriptName + "\n");
		}
		runAll.close();
	}

}
