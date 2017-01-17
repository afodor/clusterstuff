/*
 * use blast to align all samples to sRNATarBase
 */
package kw_jobinMiRNA;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;

public class BlastTosRNATarBase {
	public static String DIR = "/nobackup/afodor_research/kwinglee/jobin/microRNA/";
	public static String FQDIR = DIR + "adapterFiltered/";
	public static String OUTDIR = DIR + "sRNATarBaseBlast/";
	public static String SCRIPTDIR = "/projects/afodor_research/kwinglee/scripts/jobin/microRNA/sRNATarBase/";
	public static String REF = "/nobackup/afodor_research/kwinglee/sRNATarBase3.0/sRNATarBase3Blast";
	
	public static void main(String[] args) throws IOException {
		File odir = new File(OUTDIR);
		if(!odir.exists()) {
			odir.mkdirs();
		}
		File[] files = new File(FQDIR).listFiles();
		BufferedWriter all = new BufferedWriter(new FileWriter(new File(
				SCRIPTDIR + "blastAlignAll.sh")));
		for(File fa : files) {
			if(fa.getName().endsWith(".fasta")){
				String id = fa.getName().replace(".fasta", "");
				String scriptName = "blastAlignsRNA_" + id;
				String name = id + ".sRNAbact.blast";

				//add to run all script
				all.write("qsub -q \"copperhead\" " + scriptName + "\n");

				//write script
				BufferedWriter script = new BufferedWriter(new FileWriter(new File(
						SCRIPTDIR + scriptName)));
				script.write("#PBS -l procs=1,mem=10GB\n");
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
