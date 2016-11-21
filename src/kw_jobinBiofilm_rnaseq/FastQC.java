/*
 * run fastqc on raw reads
 */
package kw_jobinBiofilm_rnaseq;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;

public class FastQC {
	private static String FASTQC = "/nobackup/afodor_research/kwinglee/software/FastQC/fastqc";
	private static String BASEDIR = "/nobackup/afodor_research/kwinglee/jobin/biofilm/";
	private static String FQDIR = BASEDIR + "RNAseqTestRunFastqs/";
	private static String OUTDIR = BASEDIR + "rnaseq/fastqc/";
	private static String SCRIPTDIR = "/projects/afodor_research/kwinglee/scripts/jobin/biofilmRNAseq/fastqc/";
	
	public static void main(String[] args) throws IOException {
		BufferedWriter runAll = new BufferedWriter(new FileWriter(new File(
				SCRIPTDIR + "runAll.sh")));
		String[] samples = new File(FQDIR).list();
		for(String s : samples) {
			String[] fqs = new File(FQDIR + s).list();
			for(String f : fqs) {
				String name = f.replace("-L001_S1_L001", "").replace("_001.fastq.gz", "");
				String scriptName = "qc_"+ name;
				BufferedWriter script = new BufferedWriter(new FileWriter(new File(
						SCRIPTDIR + scriptName)));
				script.write(FASTQC + " --outdir=" + OUTDIR + name + " "
						+ FQDIR + s + File.separator + f + "\n");
				script.close();
				runAll.write("qsub -q \"copperhead\" " + scriptName + "\n");
				
			}
			
		}
		
		runAll.close();
	}

}
