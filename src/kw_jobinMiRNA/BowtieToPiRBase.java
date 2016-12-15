/*
 * use bowtie2 to align all samples to piRBase
 */
package kw_jobinMiRNA;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;

public class BowtieToPiRBase {
	public static String DIR = "/nobackup/afodor_research/kwinglee/jobin/microRNA/";
	public static String FQDIR = DIR + "adapterFiltered/";
	public static String OUTDIR = DIR + "piRBaseBowtie/";
	public static String SCRIPTDIR = "/projects/afodor_research/kwinglee/scripts/jobin/microRNA/piRNA/";
	public static String REF = "/nobackup/afodor_research/kwinglee/piRBase_v1.0/piRbaseMouseBowtie";
	
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
				String id = fq.getName().replace(".fastq", "");
				String scriptName = "bowtieAlignPiR_" + id;
				String name = id + ".piR.bowtie";

				//add to run all script
				all.write("qsub -q \"copperhead\" " + scriptName + "\n");

				//write script
				BufferedWriter script = new BufferedWriter(new FileWriter(new File(
						SCRIPTDIR + scriptName)));
				script.write("#PBS -l procs=1,mem=20GB,walltime=12:00:00\n");
				script.write("module load bowtie2\n");
				script.write("module load samtools\n");
				//align
				script.write("bowtie2 -x " + REF + " -U " + fq.getAbsolutePath()
						+ " -S " + OUTDIR + name + ".sam\n");
				script.write("samtools flagstat " + OUTDIR + name + ".sam > " 
						+ OUTDIR + name + ".flagstat\n");
				script.close();				
			}
		}
		all.close();
	}
}
