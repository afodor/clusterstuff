/*
 * use blast to align all samples to sRNATarBase
 */
package kw_jobinMiRNA;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;

public class BlastToSilva {
	public static String DIR = "/nobackup/afodor_research/kwinglee/jobin/microRNA/";
	public static String FQDIR = DIR + "adapterFiltered/";
	public static String OUTDIR = DIR + "silvaBlast/";
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
				SCRIPTDIR + "blastAlignAll.sh")));
		for(File fa : files) {
			if(fa.getName().endsWith(".fasta")){
				String id = fa.getName().replace(".fasta", "");
				String scriptName = "blastAlignSilva_" + id;
				String name = id + ".silva.bowtie";

				//add to run all script
				all.write("qsub -q \"copperhead\" " + scriptName + "_LSU\n");
				all.write("qsub -q \"copperhead\" " + scriptName + "_SSU\n");

				//write script LSU
				BufferedWriter script = new BufferedWriter(new FileWriter(new File(
						SCRIPTDIR + scriptName + "_LSU")));
				script.write("#PBS -l procs=1,mem=10GB\n");
				script.write("module load blast/2.5.0+\n");
				//align 
				script.write("blastn -outfmt 6 -db " + LSU + 
						" -query " + fa.getAbsolutePath()
						+ " -out " + OUTDIR + name + ".lsu.txt\n");
				script.close();
				
				//write script SSU
				script = new BufferedWriter(new FileWriter(new File(
						SCRIPTDIR + scriptName + "_SSU")));
				script.write("#PBS -l procs=1,mem=10GB\n");
				script.write("module load blast/2.5.0+\n");
				//align 
				script.write("blastn -outfmt 6 -db " + SSU + 
						" -query " + fa.getAbsolutePath()
						+ " -out " + OUTDIR + name + ".ssu.txt\n");
				script.close();
			}
		}

		all.close();
	}
}
