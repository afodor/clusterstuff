/*
 * Use HISAT2 to align all reads to mouse genome, using paired end info
 */

package kw_jobinBiofilm_rnaseq;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;

public class HisatAlignToMouse {
	public static String FQDIR = "/nobackup/afodor_research/kwinglee/jobin/biofilm/RNAseqTestRunFastqs/";
	public static String SCRIPTDIR = "/projects/afodor_research/kwinglee/scripts/jobin/biofilmRNAseq/";
	public static String OUTDIR = "/nobackup/afodor_research/kwinglee/jobin/biofilm/rnaseq/hisatAlignToMouse/";
	public static String HISATDIR = "/nobackup/afodor_research/kwinglee/software/hisat2-2.0.4/";
	public static String REF = "/nobackup/afodor_research/kwinglee/mouseGRCm38/mouseGRCm38";
	public static int NUM_THREADS = 2;

	public static void main(String[] args) throws IOException {
		BufferedWriter runAll = new BufferedWriter(new FileWriter(new File(
				SCRIPTDIR + "runHisatAlignToMouse.sh")));
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
					SCRIPTDIR + "hisatAlignToMouse_" + name)));
			//load needed modules
			script.write("#PBS -l walltime=100:00:00\n");
			script.write("#PBS -l mem=10GB\n");
			script.write("#PBS -l procs=" + NUM_THREADS + "\n");

			//align
			script.write(HISATDIR + "hisat2 -q -p " + NUM_THREADS + 
					" -x " + REF + " -1 " + FQDIR + d + File.separator + fqs[r1] + 
					" -2 " + FQDIR + d + File.separator + fqs[r2] +
					" -S " + OUTDIR + name +".hisatMouse.sam\n");

			//add to run all
			runAll.write("qsub -q \"copperhead\" hisatAlignToMouse_" + name + "\n");

			script.close();
		}
		runAll.close();
	}

}
