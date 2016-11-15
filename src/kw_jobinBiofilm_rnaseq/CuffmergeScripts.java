/*
 * make cuffmerge input tables and scripts
 */
package kw_jobinBiofilm_rnaseq;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;

public class CuffmergeScripts {
	public static String SCRIPTDIR = "/projects/afodor_research/kwinglee/scripts/jobin/biofilmRNAseq/cuffScripts/";
	public static String CUFFDIR = "/nobackup/afodor_research/kwinglee/jobin/biofilm/rnaseq/cufflinksMouseResults/";
	//public static String CUFFMERGE = "/nobackup/afodor_research/kwinglee/software/cufflinks-2.2.1.Linux_x86_64/cuffmerge";
	public static String MOUSEDIR = "/nobackup/afodor_research/kwinglee/mouseGRCm38/";
	public static String MOUSEGFF = MOUSEDIR + "GCF_000001635.25_GRCm38.p5_genomic.gff";
	public static String MOUSEFA = MOUSEDIR + "GCF_000001635.25_GRCm38.p5_genomic.fna";
	
	public static void main(String[] args) throws IOException {
		//set up manifest files
		String hisatList = CUFFDIR + "hisat_assembly_list.txt";
		BufferedWriter hisat = new BufferedWriter(new FileWriter(new File(
				hisatList)));
		String tophatList = CUFFDIR + "tophat_assembly_list.txt";
		BufferedWriter tophat = new BufferedWriter(new FileWriter(new File(
				tophatList)));
		//make manifest files
		File[] files = new File(CUFFDIR).listFiles();
		for(File f : files) {
			if(f.getName().startsWith("hisat_RNA") && f.isDirectory()) {
				hisat.write(f.getAbsolutePath() + File.separator + "transcripts.gtf\n");
			} else if(f.getName().startsWith("tophat_RNA") && f.isDirectory()) {
				tophat.write(f.getAbsolutePath() + File.separator + "transcripts.gtf\n");
			}
		}
		hisat.close();
		tophat.close();
		
		//set up scripts (for Cobra)
		makeScript("cuffmerge_tophat", tophatList);
		makeScript("cuffmerge_hisat", hisatList);
	}

	//makes the script to run cuffmerge on the given gtfList; script has scriptName as its name
	//output is to a file with scriptName as its prefix
	private static void makeScript(String scriptName, String gtfList) throws IOException {
		BufferedWriter script = new BufferedWriter(new FileWriter(new File(
				SCRIPTDIR + scriptName)));
		script.write("PATH=$PATH:/nobackup/afodor_research/kwinglee/software/cufflinks-2.2.1.Linux_x86_64\n");
		script.write("cuffmerge -o " + CUFFDIR + scriptName + " -g " + MOUSEGFF
				+ " -p 2 -s " + MOUSEFA + " " + gtfList + "\n");
		script.close();
	}
}
