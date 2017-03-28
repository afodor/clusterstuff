/**
 * generate mapping file for QIIME (run1)
 */
package kw_jobinGA;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Arrays;

public class MakeQiimeMap {
	public static String DIR = "/nobackup/afodor_research/kwinglee/jobin/ga-stool/";
	public static String FASTA_FOLDER = DIR + "sequences/";
	
	public static void main(String[] args) throws IOException {
		String[] files = new File(FASTA_FOLDER).list();
		Arrays.sort(files);
		BufferedWriter out = new BufferedWriter(new FileWriter(new File(DIR + "qiime/qiimeMap.txt")));
		out.write("#SampleID\tBarcodeSequence\tLinkerPrimerSequence\tInputFileName\tDescription\n");
		for(String name : files) {
			if(name.contains("_R1_") && name.endsWith(".fasta")) { 
				String id = name.split("_")[2].replace(".fasta", "");
				if(!id.startsWith("S")) {//include controls and gastric aspirate but not stool
					out.write(id + "\tX\tX\t" + name + "\n");					
				}
			}
		}
		out.close();
	}
}
