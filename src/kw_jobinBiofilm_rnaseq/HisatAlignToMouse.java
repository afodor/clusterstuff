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
			String name = fqs[0].replace("_001.fastq.gz", "");
			BufferedWriter script = new BufferedWriter(new FileWriter(new File(
					SCRIPTDIR + "hsatAlignToMouse_" + name)));
			//load needed modules
			script.write("#PBS -l walltime=100:00:00\n");
			script.write("#PBS -l mem=10GB\n");
			script.write("#PBS -l procs=" + NUM_THREADS + "\n");

			//align
			script.write(HISATDIR + "hisat2 -q -p " + NUM_THREADS + 
					" -x " + REF + " -1 " + fqs[0] + " -2 " + fqs[1] +
					" -S " + OUTDIR + name +".sam"
					+ "\n");

			//add to run all
			runAll.write("qsub -q \"copperhead\" hsatAlignToMouse_" + name + "\n");

			script.close();
		}
		runAll.close();
	}

}
