/*
 * Use HISAT2 to align all reads to mouse genome, using paired end info
 */

package kw_jobinBiofilm_rnaseq;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;

public class TophatAlignToMouse {
	public static String FQDIR = "/nobackup/afodor_research/kwinglee/jobin/biofilm/RNAseqTestRunFastqs/";
	public static String SCRIPTDIR = "/projects/afodor_research/kwinglee/scripts/jobin/biofilmRNAseq/";
	public static String OUTDIR = "/nobackup/afodor_research/kwinglee/jobin/biofilm/rnaseq/tophatAlignToMouse/";
	public static String TOPHAT = "nobackup/afodor_research/kwinglee/software/tophat-2.1.1.Linux_x86_64/tophat2";
	public static String REF = "/nobackup/afodor_research/kwinglee/mouseGRCm38/mouseGRCm38bowtie";
	public static int NUM_THREADS = 2;

	public static void main(String[] args) throws IOException {
		BufferedWriter runAll = new BufferedWriter(new FileWriter(new File(
				SCRIPTDIR + "runTophatAlignToMouse.sh")));
		String[] dirs = new File(FQDIR).list();
		for(String d : dirs) {
			String[] fqs = new File(FQDIR + d).list();
			int r1 = 0;//index of read1
			int r2 = 1;//index of read2
			if(!fqs[r1].contains("_R1_")) {
				r1 = 1;
				r2 = 0;
			}
			String name = fqs[r1].replaceFirst("_S[1-6]_L001_R1_001.fastq.gz", "");
			BufferedWriter script = new BufferedWriter(new FileWriter(new File(
					SCRIPTDIR + "tophatAlignToMouse_" + name)));
			//copperhead setup
			/*script.write("#PBS -l walltime=50:00:00\n");
			script.write("#PBS -l mem=10GB\n");
			script.write("#PBS -l procs=" + NUM_THREADS + "\n");*/
			
			//load needed modules
			script.write("	module load bowtie2\n");

			//align
			script.write(TOPHAT + " -o " + OUTDIR + " -p " + NUM_THREADS +
					" " + FQDIR + d + File.separator + fqs[r1] + 
					" " + FQDIR + d + File.separator + fqs[r2] + "\n");

			//add to run all
			runAll.write("qsub -q \"Cobra\" tophatAlignToMouse_" + name + "\n");

			script.close();
		}
		runAll.close();
	}

}
