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
	
	public static void main(String[] args) throws IOException {
		BufferedWriter out = new BufferedWriter(new FileWriter(new File(
				SCRIPT_DIR + "unzipHMP")));
		out.write("#PBS -l procs=1\n");
		out.write("cd " + FASTQ_DIR + "\n");
		String[] files = new File(FASTQ_DIR).list();
		for(String f : files) {
			if(f.endsWith(".tar.bz2")) {
				out.write("tar -vxjf " + f + "\n");
			}
		}
		out.close();
	}

}
