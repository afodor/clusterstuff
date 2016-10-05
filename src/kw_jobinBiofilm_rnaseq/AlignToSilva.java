/*
 * align reads to mouse genome
 */
package kw_jobinBiofilm_rnaseq;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;

public class AlignToSilva {
	public static String SSU = "/nobackup/afodor_research/kwinglee/software/silva/SILVA_128_SSURef_Nr99_tax_silva.fasta";//small subunit
	public static String LSU = "/nobackup/afodor_research/kwinglee/software/silva/SILVA_128_LSURef_tax_silva.fasta";//large subunit
	public static String FQDIR = "/nobackup/afodor_research/kwinglee/jobin/biofilm/rnaseq/mouseFilteredFastq/";
	public static String SCRIPTDIR = "/projects/afodor_research/kwinglee/scripts/jobin/biofilmRNAseq/";
	public static String OUTDIR = "/nobackup/afodor_research/kwinglee/jobin/biofilm/rnaseq/alignToSilva/";

	public static void main(String[] args) throws IOException {
		BufferedWriter runAll = new BufferedWriter(new FileWriter(new File(
				SCRIPTDIR + "runAlignToSilva.sh")));
		String[] fqs = new File(FQDIR).list();
		for(String f : fqs) {
			String name = f.replace(".fastq", "");
			BufferedWriter script = new BufferedWriter(new FileWriter(new File(
					SCRIPTDIR + "alignToSilva_" + name)));
			//load needed modules
			script.write("#PBS -l walltime=100:00:00\n");
			script.write("#PBS -l mem=6GB\n");
			script.write("module load bwa\n");
			script.write("module load samtools\n");

			//small subunit
			//align
			script.write("bwa mem " + SSU + 
					" " + FQDIR +  f + " > " +
					OUTDIR + "ssu.aln." + name + ".sam\n");//command to align file
			//get mapped reads only
			script.write("samtools view -h -S -F 4 " + OUTDIR + "ssu.aln." + name + ".sam > " 
					+ OUTDIR + "ssu.aln." + name + ".mapped.sam\n"); //get mapped reads; use -f 4 for unmapped
			//get stats
			script.write("samtools flagstat " + OUTDIR + "ssu.aln." + name + ".sam > "
					+ OUTDIR + "ssu.aln." + name + ".flagstat\n");

			//large subunit
			//align
			script.write("bwa mem " + LSU + 
					" " + FQDIR +  f + " > " +
					OUTDIR + "lsu.aln." + name + ".sam\n");//command to align file
			//get mapped reads only
			script.write("samtools view -h -S -F 4 " + OUTDIR + "lsu.aln." + name + ".sam > " 
					+ OUTDIR + "lsu.aln." + name + ".mapped.sam\n"); //get mapped reads; use -f 4 for unmapped
			//get stats
			script.write("samtools flagstat " + OUTDIR + "lsu.aln." + name + ".sam > "
					+ OUTDIR + "lsu.aln." + name + ".flagstat\n");

			//add to run all
			runAll.write("qsub -q \"copperhead\" alignToSilva_" + name + "\n");

			script.close();
		}
		runAll.close();
	}

}
