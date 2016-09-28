/*
 * generate script to combine mouse chromosomes and index them for bwa
 */

package kw_jobinBiofilm_rnaseq;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Arrays;
import java.util.Collections;

public class IndexMM10 {
	public static void main(String[] args) throws IOException {
		String mmdir = "/nobackup/afodor_research/kwinglee/mm10/";
		BufferedWriter script = new BufferedWriter(new FileWriter(new File(
				"/projects/afodor_research/kwinglee/scripts/jobin/biofilmRNAseq/indexMM10")));
		script.write("module load bwa\n");
		script.write("cd " + mmdir + "\n");
		String[] files = new File(mmdir).list();
		Arrays.sort(files);
		script.write("cat");
		for(String f : files) {
			if(f.endsWith(".fa") && !f.contains("_") && !f.contains("random")) {
				script.write(" " + f);
			}
		}
		script.write(" > mm10.fa\n");
		script.write("bwa index mm10.fa\n");
		script.close();
	}

}
