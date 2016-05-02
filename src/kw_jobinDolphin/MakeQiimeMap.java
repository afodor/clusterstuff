/**
 * generate mapping file for QIIME
 */
package kw_jobinDolphin;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;

public class MakeQiimeMap {
	public static String DIR = "/nobackup/afodor_research/kwinglee/jobin/dolphin/";
	
	public static void main(String[] args) throws IOException {
		File[] files = new File(DIR + "stitched_reads").listFiles();
		BufferedWriter out = new BufferedWriter(new FileWriter(new File(DIR + "qiimeMap.txt")));
		out.write("#SampleID\tBarcodeSequence\tLinkerPrimerSequence\tInputFileName\tDescription\n");
		for(File f : files) {
			String name = f.getName();
			if(name.endsWith("join.fasta")) {
				String id = name.replaceAll(".join.fasta", "");
				out.write(id + "\tX\tX\t" + name + "\n");
			}
		}
		out.close();
	}
}
