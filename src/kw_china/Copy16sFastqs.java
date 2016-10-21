/*
 * generate script to copy 16S fastqs, remove the .filepart from the name
 * and giving them a more meaningful name so can upload to NCBI SRA,
 * and produce list of final names
 */

package kw_china;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Arrays;

public class Copy16sFastqs {
	public static final String SEQDIR = "/nobackup/afodor_research/ChinaSequences/";
	public static final String NEWDIR = "/nobackup/afodor_research/kwinglee/china/fastqs_16s/";
	
	public static void main(String[] args) throws IOException {
		BufferedWriter script = new BufferedWriter(new FileWriter(new File(
				NEWDIR + "copyFastqs.sh")));
		BufferedWriter fileList = new BufferedWriter(new FileWriter(new File(
				NEWDIR + "fastqList.txt")));
		String[] dirs = new String[]{"first", "second"};
		for(String time : dirs) {
			String path = SEQDIR + time + "/microbiome/";
			String[] subdirs = new File(path).list();
			for(String sub : subdirs) {
				path = SEQDIR + time + "/microbiome/" + sub + "/Clean/";
				String[] sampleDirs = new File(path).list();
				for(String sample : sampleDirs) {
					String[] fqs = new File(path + sample).list();
					Arrays.sort(fqs);
					//check read 1 is first in array
					if(!fqs[0].contains("_1.fq.gz") && !fqs[1].contains("_2.fq.gz")) {
						System.err.println("Bad sort: " + path + sample);
						System.err.println("\t" + fqs[0] + "\t\n" + fqs[1]);
					}
					String newName = time + "_" + sample;
					script.write("cp " + path + sample + File.separator + fqs[0] +
							" " + newName + "_1.fq.gz\n");
					script.write("cp " + path + sample + File.separator + fqs[1] +
							" " + newName + "_2.fq.gz\n");
					fileList.write(newName + "_1.fq.gz\n" + newName + "_2.fq.gz\n");
				}
			}
		}
		script.close();
		fileList.close();
	}
}
