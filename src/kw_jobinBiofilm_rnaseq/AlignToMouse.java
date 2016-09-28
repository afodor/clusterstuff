/*
 * align reads to mouse genome
 */
package kw_jobinBiofilm_rnaseq;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;

public class AlignToMouse {
	public static String MM10 = "/nobackup/afodor_research/kwinglee/mm10/mm10.fa";
	public static String FQDIR = "/nobackup/afodor_research/kwinglee/jobin/biofilm/RNAseqTestRunFastqs/";
	public static String SCRIPTDIR = "/projects/afodor_research/kwinglee/scripts/jobin/biofilmRNAseq/";
	public static String OUTDIR = "/nobackup/afodor_research/kwinglee/jobin/biofilm/rnaseq/alignToMouse/";
	
	public static void main(String[] args) throws IOException {
		BufferedWriter runAll = new BufferedWriter(new FileWriter(new File(
				SCRIPTDIR + "runAlignToMouse.sh")));
		String[] dirs = new File(FQDIR).list();
		for(String d : dirs) {
			String[] fqs = new File(FQDIR + d).list();
			for(String f :fqs) {
				String name = f.replace("_001.fastq.gz", "");
				BufferedWriter script = new BufferedWriter(new FileWriter(new File(
						SCRIPTDIR + "alignToMouse_" + name)));
				//load needed modules
				script.write("#PBS -l walltime=100:00:00\n");
				script.write("module load bwa\n");
				script.write("module load samtools\n");
				
				//align
				script.write("bwa mem " + MM10 + 
						" " + FQDIR + d + "/" + f + " > " +
						OUTDIR + "mm10aln." + name + ".sam\n");//command to align file
				//get mapped reads only
				script.write("samtools view -h -S -F 4 " + OUTDIR + "mm10aln." + name + ".sam > " 
						+ OUTDIR + "mm10aln." + name + ".mapped.sam\n"); //get mapped reads; use -f 4 for unmapped
				
				//add to run all
				runAll.write("qsub -q \"copperhead\" alignToMouse_" + name + "\n");
				
				script.close();
			}
		}
		runAll.close();
	}

}
