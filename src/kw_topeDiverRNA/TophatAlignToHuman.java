/*
 * Use Tophat2 to align all reads to human genome, using paired end info
 */

package kw_topeDiverRNA;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;

public class TophatAlignToHuman {
	public static String FQDIR = "/nobackup/afodor_research/topeRnaSeqFeb2017/170118_UNC31-K00269_0051_AHGWY3BBXX/";
	public static String SCRIPTDIR = "/projects/afodor_research/kwinglee/scripts/tope/diverRNA/";
	public static String OUTDIR = "/nobackup/afodor_research/kwinglee/tope/diverRNAseq/tophatAlignToHuman/";
	public static String TOPHAT = "/nobackup/afodor_research/kwinglee/software/tophat-2.1.1.Linux_x86_64/tophat2";
	public static String REF = "/nobackup/afodor_research/kwinglee/hg38/hg38bowtie";
	public static int NUM_THREADS = 2;

	public static void main(String[] args) throws IOException {
		BufferedWriter runAll = new BufferedWriter(new FileWriter(new File(
				SCRIPTDIR + "runTophatAlignToHuman.sh")));
		String[] fqs = new File(FQDIR).list();
		for(String fq : fqs) {
			if(fq.endsWith("_R1_001.fastq.gz")) {
				String name = fq.split("_")[0].replaceAll("-$", "");
				String scriptName = "tophatHuman_" + name;
				BufferedWriter script = new BufferedWriter(new FileWriter(new File(
						SCRIPTDIR + scriptName)));
				//copperhead setup
				script.write("#PBS -l walltime=600:00:00\n");
				/*script.write("#PBS -l mem=10GB\n");
				script.write("#PBS -l procs=" + NUM_THREADS + "\n");*/

				//load needed modules
				script.write("module load bowtie2\n");

				//align
				script.write(TOPHAT + " -o " + OUTDIR + name + " -p " + NUM_THREADS +
						" " + REF + " " + FQDIR + fq + 
						" " + FQDIR + fq.replace("_R1_", "_R2_") + "\n");

				//add to run all
				runAll.write("qsub -q \"Cobra\" " + scriptName + "\n");

				script.close();
			}
		}
		runAll.close();
	}

}
