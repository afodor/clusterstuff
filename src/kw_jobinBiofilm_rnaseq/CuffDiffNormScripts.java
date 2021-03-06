/*
 * generate scripts to run cuffdiff and cuffnorm
 */
package kw_jobinBiofilm_rnaseq;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;

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
		String indiv = "";
		String indLab = "";
		for(File f : files) {
			String name = f.getName();
			if(name.contains("M115") || name.contains("M317") || name.contains("M426")) {//il10
				il10 += f.getAbsolutePath() + File.separator + "accepted_hits.bam,";
			} else if(name.contains("M527") || name.contains("M635") || name.contains("M937")) { //apcmin
				apc += f.getAbsolutePath() + File.separator + "accepted_hits.bam,";
			} else {
				throw new Exception("Invalid name: " + name);
			}
			indiv += f.getAbsolutePath() + File.separator + "accepted_hits.bam ";
			indLab += name + ",";
		}
		il10 = il10.replaceAll(",$", "");
		apc = apc.replaceAll(",$", "");
		indiv = indiv.replaceAll(" $", "");
		indLab = indLab.replaceAll(",$", "");
		
		//cuffdiff
		BufferedWriter diff = new BufferedWriter(new FileWriter(new File(
				SCRIPTDIR + "cuffdiff_tophat")));
		diff.write("PATH=$PATH:" + CUFFDIR + "\n");
		diff.write("cuffdiff -o " + OUTDIR + "cuffdiff_tophat" 
				+ " -L ApcMinIL10KO,ApcMin -p 2 " + GFFMERGE + " " +
				il10 + " " + apc + "\n");
		diff.close();
		
		//cuffnorm by genotype
		BufferedWriter norm = new BufferedWriter(new FileWriter(new File(
				SCRIPTDIR + "cuffnorm_tophat")));
		norm.write("PATH=$PATH:" + CUFFDIR + "\n");
		norm.write("cuffnorm -o " + OUTDIR + "cuffnorm_tophat" 
				+ " -L ApcMinIL10KO,ApcMin -p 2 -library-norm-method classic-fpkm " + 
				GFFMERGE + " " + il10 + " " + apc + "\n");
		norm.close();
		
		//cuffnorm separated
		BufferedWriter norm2 = new BufferedWriter(new FileWriter(new File(
				SCRIPTDIR + "cuffnorm_tophat_indiv")));
		norm2.write("PATH=$PATH:" + CUFFDIR + "\n");
		norm2.write("cuffnorm -o " + OUTDIR + "cuffnorm_tophat_indiv" 
				+ " -L " + indLab + " -p 4 -library-norm-method classic-fpkm " + 
				GFFMERGE + " " + indiv + "\n");
		norm2.close();
		
		//cuffdiff for cages
		String cage1 = "";
		String cage2 = "";
		String cage3 = "";
		for(File f : files) {
			String name = f.getName();
			if(name.contains("M115") || name.contains("M527") || name.contains("M635")) {//cage1
				cage1 += f.getAbsolutePath() + File.separator + "accepted_hits.bam,";
			} else if(name.contains("M937")) { //cage2
				cage2 += f.getAbsolutePath() + File.separator + "accepted_hits.bam,";
			} else if(name.contains("M317") || name.contains("M426")) { //cage3
				cage3 += f.getAbsolutePath() + File.separator + "accepted_hits.bam,";
			} else {
				throw new Exception("Invalid name: " + name);
			}
		}
		cage1 = cage1.replaceAll(",$", "");
		cage2 = cage2.replaceAll(",$", "");
		cage3 = cage3.replaceAll(",$", "");
		BufferedWriter cage = new BufferedWriter(new FileWriter(new File(
				SCRIPTDIR + "cuffdiff_tophat_cage")));
		cage.write("PATH=$PATH:" + CUFFDIR + "\n");
		cage.write("cuffdiff -o " + OUTDIR + "cuffdiff_tophat_cage" 
				+ " -L cage1,cage2,cage3 -p 2 " + GFFMERGE + " " +
				cage1 + " " + cage2 + " " + cage3 + "\n");
		cage.close();
	}

}
