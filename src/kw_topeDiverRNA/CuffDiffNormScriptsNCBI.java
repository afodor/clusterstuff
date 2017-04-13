/*
 * generate scripts to run cuffdiff and cuffnorm
 */
package kw_topeDiverRNA;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.util.HashMap;

public class CuffDiffNormScriptsNCBI {
	private static String CUFFDIR = "/nobackup/afodor_research/kwinglee/software/cufflinks-2.2.1.Linux_x86_64";
	private static String SCRIPTDIR = "/projects/afodor_research/kwinglee/scripts/tope/diverRNA/cufflinks/";
	private static String BASEDIR = "/nobackup/afodor_research/kwinglee/tope/diverRNAseq/";
	private static String OUTDIR = BASEDIR + "cufflinksHumanResultsNCBI/";
	private static String GFFMERGE = BASEDIR + "cufflinksHumanResultsNCBI/cuffmergeNCBI/merged.gtf";

	public static void main(String[] args) throws Exception {
		//get CaseControl
		BufferedReader meta = new BufferedReader(new FileReader(new File(BASEDIR + "RNAseqMetadata.txt")));
		HashMap<String, Integer> caseControl = new HashMap<String, Integer>();
		int ccCol = 14;
		int idCol = 0;
		//get caseControl column
		String header = meta.readLine();
		String[] sp = header.split("\t");
		if(!sp[idCol].equals("StudyID")) {
			meta.close();
			throw new Exception("Incorrect ID column");
		}
		if(!sp[ccCol].equals("CaseControl")) {
			meta.close();
			throw new Exception("Incorrect case control column");
		}
		for(String line = meta.readLine(); line != null; line = meta.readLine()) {
			sp = line.split("\t");
			if(sp.length < ccCol) {
				caseControl.put(sp[idCol], 2);//control missing annotation
			} else {
				caseControl.put(sp[idCol], Integer.parseInt(sp[ccCol]));	
			}
		}
		meta.close();

		//get list of aligned files
		File[] files = new File(BASEDIR + "tophatAlignToHumanNCBI/").listFiles();
		String cc0 = "";
		String cc1 = "";
		String indiv = "";
		String indLab = "";
		for(File f : files) {
			String name = f.getName();
			if(!caseControl.containsKey(name)) {
				throw new Exception("Missing sample: " + name);
			} else if(caseControl.get(name) == 0) { 
				cc0 += f.getAbsolutePath() + File.separator + "accepted_hits.bam,";
			} else if(caseControl.get(name) == 1) {
				cc1 += f.getAbsolutePath() + File.separator + "accepted_hits.bam,";
			} else {
				System.err.println("Skipped caseControl: " + name + " " + caseControl.get(name));
			}
			indiv += f.getAbsolutePath() + File.separator + "accepted_hits.bam ";
			indLab += name + ",";
		}
		cc0 = cc0.replaceAll(",$", "");
		cc1 = cc1.replaceAll(",$", "");
		indiv = indiv.replaceAll(" $", "");
		indLab = indLab.replaceAll(",$", "");

		//cuffdiff
		BufferedWriter diff = new BufferedWriter(new FileWriter(new File(
				SCRIPTDIR + "cuffdiffNCBI")));
		diff.write("#PBS -l procs=1,walltime=150:00:00\n");
		diff.write("PATH=$PATH:" + CUFFDIR + "\n");
		diff.write("cuffdiff -o " + OUTDIR + "cuffdiffNCBI_caseControl" 
				+ " -L 0,1 -p 2 " + GFFMERGE + " " +
				cc0 + " " + cc1 + "\n");
		diff.close();

		//cuffnorm by case control
		BufferedWriter norm = new BufferedWriter(new FileWriter(new File(
				SCRIPTDIR + "cuffnormNCBI_caseControl")));
		norm.write("#PBS -l procs=1,walltime=150:00:00\n");
		norm.write("PATH=$PATH:" + CUFFDIR + "\n");
		norm.write("cuffnorm -o " + OUTDIR + "cuffnormNCBI_caseControl" 
				+ " -L 0,1 -p 2 -library-norm-method classic-fpkm " + 
				GFFMERGE + " " + cc0 + " " + cc1 + "\n");
		norm.close();

		//cuffnorm separated
		BufferedWriter norm2 = new BufferedWriter(new FileWriter(new File(
				SCRIPTDIR + "cuffnormNCBI_indiv")));
		norm2.write("#PBS -l procs=1,walltime=150:00:00\n");
		norm2.write("PATH=$PATH:" + CUFFDIR + "\n");
		norm2.write("cuffnorm -o " + OUTDIR + "cuffnormNCBI_indiv" 
				+ " -L " + indLab + " -p 4 -library-norm-method classic-fpkm " + 
				GFFMERGE + " " + indiv + "\n");
		norm2.close();
	}

}
