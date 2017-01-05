/*
 * use bowtie2 to align all samples to Silva large and small subunit
 */
package kw_jobinMiRNA;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;

public class BowtieToSilva {
	public static String DIR = "/nobackup/afodor_research/kwinglee/jobin/microRNA/";
	public static String FQDIR = DIR + "adapterFiltered/";
	public static String OUTDIR = DIR + "silvaBowtie/";
	public static String SCRIPTDIR = "/projects/afodor_research/kwinglee/scripts/jobin/microRNA/silva/";
	public static String SSU = "/nobackup/afodor_research/kwinglee/software/silva/SILVA_128_SSURef_Nr99_tax_silva";
	public static String LSU = "/nobackup/afodor_research/kwinglee/software/silva/SILVA_128_LSURef_tax_silva";
	
	public static void main(String[] args) throws IOException {
		File odir = new File(OUTDIR);
		if(!odir.exists()) {
			odir.mkdirs();
		}
		File[] files = new File(FQDIR).listFiles();
		BufferedWriter all = new BufferedWriter(new FileWriter(new File(
				SCRIPTDIR + "bowtieAlignAllSilva.sh")));
		for(File fq : files) {
			if(fq.getName().endsWith(".fastq")) {
				String id = fq.getName().replace(".fastq", "");
				String scriptName = "bowtieAlignSilva_" + id;
				String name = id + ".silva.bowtie";

				//add to run all script
				all.write("qsub -q \"copperhead\" " + scriptName + "_LSU\n");
				all.write("qsub -q \"copperhead\" " + scriptName + "_SSU\n");

				//write script for LSU
				BufferedWriter script = new BufferedWriter(new FileWriter(new File(
						SCRIPTDIR + scriptName + "_LSU")));
				script.write("#PBS -l procs=1,mem=10GB\n");
				script.write("module load bowtie2\n");
				script.write("module load samtools\n");
				//align 
				script.write("bowtie2 -x " + LSU + " -U " + fq.getAbsolutePath()
						+ " -S " + OUTDIR + name + ".lsu.sam\n");
				script.write("samtools flagstat " + OUTDIR + name + ".lsu.sam > " 
						+ OUTDIR + name + ".lsu.flagstat\n");
				
				//write script for SSU
				script = new BufferedWriter(new FileWriter(new File(
						SCRIPTDIR + scriptName + "_SSU")));
				script.write("#PBS -l procs=1,mem=10GB\n");
				script.write("module load bowtie2\n");
				script.write("module load samtools\n");
				//align 
				script.write("bowtie2 -x " + SSU + " -U " + fq.getAbsolutePath()
						+ " -S " + OUTDIR + name + ".ssu.sam\n");
				script.write("samtools flagstat " + OUTDIR + name + ".ssu.sam > " 
						+ OUTDIR + name + ".ssu.flagstat\n");
				script.close();
			}
		}
		all.close();
	}
}
