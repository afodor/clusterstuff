/*
 * generate the scripts to run cufflinks on Tophat results
 */
package kw_jobinBiofilm_rnaseq;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;

public class CufflinksTophatScripts {
	private static String BASEDIR = "/nobackup/afodor_research/kwinglee/jobin/biofilm/rnaseq/";
	private static String TOPHATDIR = BASEDIR + "tophatAlignToMouse/";
	private static String CUFF = "/nobackup/afodor_research/kwinglee/software/cufflinks-2.2.1.Linux_x86_64/cufflinks";
	private static String OUTDIR = BASEDIR + "cufflinksMouseResults/";
	private static String SCRIPTDIR = "/projects/afodor_research/kwinglee/scripts/jobin/biofilmRNAseq/cuffScripts/";
	private static String GFF = "/nobackup/afodor_research/kwinglee/mouseGRCm38/GCF_000001635.25_GRCm38.p5_genomic.gff";

	public static void main(String[] args) throws IOException {
		String[] sams = new File(TOPHATDIR).list();
		BufferedWriter runAll = new BufferedWriter(new FileWriter(new File(
				SCRIPTDIR + "runCufflinksTophat.sh")));
		for(String s : sams) {
			String scriptName = "cufflinksTophat_" + s; 
			BufferedWriter script = new BufferedWriter(new FileWriter(new File(
					SCRIPTDIR + scriptName)));
			/*script.write("#PBS -l mem=4GB\n");
			script.write("#PBS -l walltime=50:00:00\n");*/
			script.write(CUFF + " --GTF-guide " + GFF + " --output-dir " + 
					OUTDIR + "tophat_" + s + " " + TOPHATDIR + s + 
					File.separator + "accepted_hits.bam\n");
			script.close();

			//runAll.write("qsub -q \"copperhead\" " + scriptName + "\n");
			runAll.write("qsub -q \"Cobra\" " + scriptName + "\n");
		}
		runAll.close();
	}

}
