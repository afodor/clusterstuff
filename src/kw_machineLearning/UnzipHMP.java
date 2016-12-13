/*
 * Need to unzip HMP files before running kraken
 */
package kw_machineLearning;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;

public class UnzipHMP {
	public static String FASTQ_DIR = "/nobackup/afodor_research/kwinglee/machineLearning/hmp/fastqs/stool/";
	public static String SCRIPT_DIR = "/projects/afodor_research/kwinglee/scripts/machineLearning/hmp/";
	public static int MAXCMDS = 5;
	
	public static void main(String[] args) throws IOException {
		int numCmds = 0;
		int numScripts = 0;
		String sbase = "unzipHMP_";
		BufferedWriter all = new BufferedWriter(new FileWriter(new File(
				SCRIPT_DIR + "unzipAll.sh")));
		all.write("qsub -q \"copperhead\" " + sbase + numScripts + "\n");
		BufferedWriter script = new BufferedWriter(new FileWriter(new File(
				SCRIPT_DIR + sbase + numScripts)));
		script.write("#PBS -l procs=1\n");
		script.write("cd " + FASTQ_DIR + "\n");
		String[] files = new File(FASTQ_DIR).list();
		for(String f : files) {
			if(f.endsWith(".tar.bz2")) {
				script.write("tar -vxjf " + f + "\n");
				numCmds++;
				if(numCmds == MAXCMDS) {
					numCmds = 0;
					numScripts++;
					script.close();
					script = new BufferedWriter(new FileWriter(new File(
							SCRIPT_DIR + sbase + numScripts)));
					script.write("#PBS -l procs=1\n");
					script.write("cd " + FASTQ_DIR + "\n");
					all.write("qsub -q \"copperhead\" " + sbase + numScripts + "\n");
				}
			}
		}
		script.close();
		all.close();
		System.out.println("Number of scripts: " + (numScripts+1));
	}

}
