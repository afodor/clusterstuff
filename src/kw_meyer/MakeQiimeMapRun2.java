/**
 * generate mapping file for QIIME (run1)
 */
package kw_meyer;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Arrays;

public class MakeQiimeMapRun2 {
	public static String FASTA_FOLDER = "/nobackup/afodor_research/kwinglee/meyer/run2joinedReads/";
	public static String DIR = "/nobackup/afodor_research/kwinglee/meyer/";
	
	public static void main(String[] args) throws IOException {
		String[] files = new File(FASTA_FOLDER).list();
		Arrays.sort(files);
		BufferedWriter out = new BufferedWriter(new FileWriter(new File(DIR + "run2qiime/qiimeMapRun2.txt")));
		out.write("#SampleID\tBarcodeSequence\tLinkerPrimerSequence\tInputFileName\tDescription\n");
		for(String name : files) {
			if(name.endsWith("join.fasta")) { 
				String id = name.replace("join.fasta", "").replace("-", ".").replace("_", ".");//can't have - or _ in sampleID
				out.write(id + "\tX\tX\t" + name + "\n");
			}
		}
		out.close();
	}
}
