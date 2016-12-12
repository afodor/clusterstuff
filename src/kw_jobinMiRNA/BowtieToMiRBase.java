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
	public static String FQDIR = DIR + "fastqs/";
	public static String OUTDIR = DIR + "miRBaseAlign/";
	public static String SCRIPTDIR = "/projects/afodor_research/kwinglee/scripts/jobin/microRNA/";
	public static String REF = "/nobackup/afodor_research/kwinglee/mirbase_v21/matureBowtie";
	
	public static void main(String[] args) throws IOException {
		File[] fqs = new File(FQDIR).listFiles();
		BufferedWriter all = new BufferedWriter(new FileWriter(new File(
				SCRIPTDIR + "miRalignAll.sh")));
		for(File fq : fqs) {
			String name = "miRalign_" + fq.getName().split("-")[0];
			
			//add to run all script
			all.write("qsub -q \"copperhead\" " + name + "\n");
			
			//write script
			BufferedWriter script = new BufferedWriter(new FileWriter(new File(
					SCRIPTDIR + name)));
			script.write("module load bowtie2\n");
			script.write("bowtie2 -x " + REF + " -U " + fq.getAbsolutePath()
					+ " -S " + OUTDIR + name + "\n");
			script.close();
		}
		
		all.close();
	}
}
