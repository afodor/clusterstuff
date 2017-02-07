/*
 * set up scripts to align biofilm and APC datasets to F. rodentium
 */
package kw_jobinFrodentium;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;

public class BlastScripts {
	private static String DIR = "/nobackup/afodor_research/kwinglee/jobin/F_rodentium/";
	private static String REF = DIR + "BlastDB_Faecalibaculum_rodentium_SilvaSSU_DNA";
	private static String SCRIPTDIR = DIR + "blastScripts/";
	private static String OUTDIR = DIR + "blastResults/";

	public static void main(String[] args) throws IOException {
		//make directories
		File temp = new File(SCRIPTDIR);
		if(!temp.exists()) {
			temp.mkdirs();
		}
		temp = new File(OUTDIR);
		if(!temp.exists()) {
			temp.mkdirs();
		}
		
		//run scripts for each set of data
		String[] fastaDirs = new String[]{"/nobackup/afodor_research/kwinglee/jobin/apcTumor/fastas/",
				"/nobackup/afodor_research/kwinglee/jobin/biofilm/fastas/",
				"/nobackup/afodor_research/kwinglee/jobin/gemcitabine/demultiplexedReads/"};
		String[] names = new String[]{"apc", "biofilm", "biofilmReassoc"};
		BufferedWriter runAll = new BufferedWriter(new FileWriter(new File(
				SCRIPTDIR + "runAll.sh")));
		for(int i = 0; i < fastaDirs.length; i++) {
			String fd = fastaDirs[i];
			File[] fastas = new File(fd).listFiles();
			for(File f : fastas) {
				String n = f.getName();
				if(n.endsWith(".fasta") && ! n.contains("PCR") &&
						!n.contains("NC101") && !n.contains("H20") &&
						!n.contains("water") && !n.contains("other") &&
						(!fd.contains("gemcitabine") || n.startsWith("B"))) {
					String scriptName = "Frod_" + names[i] + "_" + n.replace(".fasta", "");
					runAll.write("qsub -q \"copperhead\" " + scriptName + "\n");
					
					BufferedWriter script = new BufferedWriter(new FileWriter(new File(
							SCRIPTDIR + scriptName)));
					script.write("#PBS -l procs=1\n");
					script.write("module load blast/2.5.0+\n");
					script.write("blastn -outfmt 7 -db " + REF + 
							" -query " + f.getAbsolutePath()
							+ " -out " + OUTDIR + scriptName + ".txt\n");
					script.close();
				}
			}
		}
		runAll.close();
	}
}
