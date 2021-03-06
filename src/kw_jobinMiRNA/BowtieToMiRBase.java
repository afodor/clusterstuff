/*
 * use bowtie2 to align all samples to miRBase
 */
package kw_jobinMiRNA;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;

public class BowtieToMiRBase {
	public static String DIR = "/nobackup/afodor_research/kwinglee/jobin/microRNA/";
	public static String FQDIR = DIR + "adapterFiltered/";
	public static String OUTDIR = DIR + "miRBaseBowtie/";
	public static String SCRIPTDIR = "/projects/afodor_research/kwinglee/scripts/jobin/microRNA/";
	public static String MATREF = "/nobackup/afodor_research/kwinglee/mirbase_v21/matureDNABowtie";
	public static String PINREF = "/nobackup/afodor_research/kwinglee/mirbase_v21/hairpinDNABowtie";

	public static void main(String[] args) throws IOException {
		File odir = new File(OUTDIR);
		if(!odir.exists()) {
			odir.mkdirs();
		}
		File[] files = new File(FQDIR).listFiles();
		BufferedWriter all = new BufferedWriter(new FileWriter(new File(
				SCRIPTDIR + "bowtieAlignAll.sh")));
		for(File fq : files) {
			if(fq.getName().endsWith(".fastq")) {
				/*String scriptName = "bowtieAlign_" + fq.getName().split("-")[0];
				String name = fq.getName().split("-")[0] + ".bowtie";*/
				String id = fq.getName().replace(".fastq", "");
				String scriptName = "bowtieAlign_" + id;
				String name = id + ".bowtie";

				//add to run all script
				all.write("qsub -q \"copperhead\" " + scriptName + "\n");

				//write script
				BufferedWriter script = new BufferedWriter(new FileWriter(new File(
						SCRIPTDIR + scriptName)));
				script.write("#PBS -l procs=1\n");
				script.write("module load bowtie2\n");
				script.write("module load samtools\n");
				//align to mature
				script.write("bowtie2 -x " + MATREF + " -U " + fq.getAbsolutePath()
						+ " -S " + OUTDIR + name + ".mature.sam\n");
				script.write("samtools flagstat " + OUTDIR + name + ".mature.sam > " 
						+ OUTDIR + name + ".mature.flagstat\n");
				//align to hairpin
				script.write("bowtie2 -x " + PINREF + " -U " + fq.getAbsolutePath()
						+ " -S " + OUTDIR + name + ".hairpin.sam\n");//align to hairpin
				script.write("samtools flagstat " + OUTDIR + name + ".hairpin.sam > " 
						+ OUTDIR + name + ".hairpin.flagstat\n");
				script.close();
			}
		}
		all.close();
	}
}
