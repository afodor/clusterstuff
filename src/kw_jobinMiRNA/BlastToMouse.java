/*
 * use blast to align all samples to miRBase
 */
package kw_jobinMiRNA;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;

public class BlastToMouse {
	public static String DIR = "/nobackup/afodor_research/kwinglee/jobin/microRNA/";
	public static String FQDIR = DIR + "adapterFiltered/";
	public static String OUTDIR = DIR + "mouseBlast/";
	public static String SCRIPTDIR = "/projects/afodor_research/kwinglee/scripts/jobin/microRNA/mouse/";
	public static String REF = "/nobackup/afodor_research/kwinglee/mm10/mm10Blast";
	
	public static void main(String[] args) throws IOException {
		File odir = new File(OUTDIR);
		if(!odir.exists()) {
			odir.mkdirs();
		}
		File[] files = new File(FQDIR).listFiles();
		BufferedWriter all = new BufferedWriter(new FileWriter(new File(
				SCRIPTDIR + "blastAlignAllMouse.sh")));
		for(File fa : files) {
			if(fa.getName().endsWith(".fasta")){
				/*String scriptName = "blastAlign_" + fa.getName().split("-")[0];
				String name = fa.getName().split("-")[0] + ".blast";*/
				String id = fa.getName().replace(".fasta", "");
				String scriptName = "blastAlign_" + id;
				String name = id + ".mouse.blast";

				//add to run all script
				all.write("qsub -q \"copperhead\" " + scriptName + "\n");

				//write script
				BufferedWriter script = new BufferedWriter(new FileWriter(new File(
						SCRIPTDIR + scriptName)));
				script.write("#PBS -l procs=1\n");
				script.write("module load blast/2.5.0+\n");
				//align 
				script.write("blastn -outfmt 6 -db " + REF + 
						" -query " + fa.getAbsolutePath()
						+ " -out " + OUTDIR + name + ".txt\n");
				script.close();
			}
		}

		all.close();
	}
}
