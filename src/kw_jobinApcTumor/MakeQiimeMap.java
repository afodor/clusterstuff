/**
 * generate mapping file for QIIME (mouse samples only)
 */
package kw_jobinApcTumor;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Arrays;

public class MakeQiimeMap {
	public static String DIR = "/nobackup/afodor_research/kwinglee/jobin/apcTumor/";
	
	public static void main(String[] args) throws IOException {
		String[] files = new File(DIR + "joinedReads").list();
		Arrays.sort(files);
		BufferedWriter out = new BufferedWriter(new FileWriter(new File(DIR + "qiime/qiimeMap.txt")));
		out.write("#SampleID\tBarcodeSequence\tLinkerPrimerSequence\tInputFileName\tDescription\n");
		for(String name : files) {
			if(name.endsWith("join.fasta") && 
					!name.contains("water") && 
					!name.contains("PCR") &&
					!name.contains("-")) { //only want the mouse samples; no controls, pancreas or biopsy
				String id = name.replace("join.fasta", "");
				out.write(id + "\tX\tX\t" + name + "\n");
			}
		}
		out.close();
	}
}
