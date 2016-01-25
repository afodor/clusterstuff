/**
 * generate mapping file for QIIME
 */
package kw_jobinAnaerobe;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;

public class makeQiimeMap {
	public static String DIR = "/nobackup/afodor_research/kwinglee/jobin/anaerobe/";
	
	public static void main(String[] args) throws IOException {
		File[] fastas = new File(DIR + "join").listFiles();
		BufferedWriter out = new BufferedWriter(new FileWriter(new File(DIR + "qiimeMap.txt")));
		out.write("#SampleID\tBarcodeSequence\tLinkerPrimerSequence\tInputFileName\tDescription\n");
		for(File f : fastas) {
			String name = f.getName();
			if(name.endsWith("join.fasta")) {
				String id = name.replaceAll("join.fasta", "").replaceAll("anaerobe_", "");
				out.write(id + "\tX\tX\t" + name + "\n");
			}
		}
		out.close();
	}
}
