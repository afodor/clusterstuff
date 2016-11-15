/*
 * generate scripts to run cuffdiff and cuffnorm
 */
package kw_jobinBiofilm_rnaseq;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;

public class CuffDiffNormScripts {
	private static String CUFFDIR = "/nobackup/afodor_research/kwinglee/software/cufflinks-2.2.1.Linux_x86_64";
	private static String SCRIPTDIR = "/projects/afodor_research/kwinglee/scripts/jobin/biofilmRNAseq/cuffScripts/";
	private static String BASEDIR = "/nobackup/afodor_research/kwinglee/jobin/biofilm/rnaseq/";
	private static String OUTDIR = BASEDIR + "cufflinksMouseResults/";
	private static String GFFMERGE = BASEDIR + "cufflinksMouseResults/cuffmerge_tophat/merged.gtf";
	
	public static void main(String[] args) throws Exception {
		//get list of aligned files
		File[] files = new File(BASEDIR + "tophatAlignToMouse/").listFiles();
		String il10 = "";
		String apc = "";
		for(File f : files) {
			String name = f.getName();
			if(name.contains("M115") || name.contains("M317") || name.contains("M426")) {//il10
				il10 += f.getAbsolutePath() + File.separator + "accepted_hits.bam,";
			} else if(name.contains("M527") || name.contains("M635") || name.contains("M937")) { //apcmin
				apc += f.getAbsolutePath() + File.separator + "accepted_hits.bam,";
			} else {
				throw new Exception("Invalid name: " + name);
			}
		}
		il10.replaceAll(",$", "");
		apc.replaceAll(",$", "");
		
		//cuffdiff
		BufferedWriter diff = new BufferedWriter(new FileWriter(new File(
				SCRIPTDIR + "cuffdiff_tophat")));
		diff.write("PATH=$PATH:" + CUFFDIR + "\n");
		diff.write("cuffdiff -o " + OUTDIR + "cuffdiff_tophat" 
				+ " -L ApcMinIL10KO, ApcMin -p 2 " + GFFMERGE + " " +
				il10 + " " + apc + "\n");
		diff.close();
		
		//cuffnorm
		BufferedWriter norm = new BufferedWriter(new FileWriter(new File(
				SCRIPTDIR + "cuffnorm_tophat")));
		norm.write("PATH=$PATH:" + CUFFDIR + "\n");
		norm.write("cuffnorm -o " + OUTDIR + "cuffnorm_tophat" 
				+ " -L ApcMinIL10KO, ApcMin -p 2 –library-norm-method classic-fpkm " + 
				GFFMERGE + " " + il10 + " " + apc + "\n");
		norm.close();
	}

}
