/*
 * generate the scripts to run cufflinks on HiSat results
 */
package kw_jobinBiofilm_rnaseq;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;

public class CufflinksScripts {
	private static String BASEDIR = "/nobackup/afodor_research/kwinglee/jobin/biofilm/rnaseq/";
	private static String HISATDIR = BASEDIR + "hisatAlignToMouse/";
	private static String CUFF = "/nobackup/afodor_research/kwinglee/software/cufflinks-2.2.1.Linux_x86_64/cufflinks";
	private static String OUTDIR = BASEDIR + "cufflinksMouseResults/";
	private static String SCRIPTDIR = "/projects/afodor_research/kwinglee/scripts/jobin/biofilmRNAseq/cuffScripts/";
	private static String GFF = "/nobackup/afodor_research/kwinglee/mouseGRCm38/GCF_000001635.25_GRCm38.p5_genomic.gff";
	
	public static void main(String[] args) throws IOException {
		String[] sams = new File(HISATDIR).list();
		BufferedWriter runAll = new BufferedWriter(new FileWriter(new File(
				SCRIPTDIR + "runCufflinks.sh")));
		for(String s : sams) {
			if(s.endsWith(".sam")) {
				String scriptName = "cufflinks_" + s.replace(".hisatMouse.sam", ""); 
				BufferedWriter script = new BufferedWriter(new FileWriter(new File(
						SCRIPTDIR + scriptName)));
				script.write("#PBS -l mem=4GB\n");
				//script.write("module load samtools\n");
				//script.write("samtools sort -o " + HISATDIR + s + ".sorted " + HISATDIR + s + "\n");
				script.write(CUFF + " --GTF-guide " + GFF + " --output-dir " + OUTDIR + " " + HISATDIR + s + ".sorted\n");
				script.close();
				
				runAll.write("qsub -q \"copperhead\" " + scriptName + "\n");
			}
		}
		runAll.close();
	}

}
